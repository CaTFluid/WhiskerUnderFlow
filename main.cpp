// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/equation_systems.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/gmv_io.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibamr/IBFECentroidPostProcessor.h>
// Elasticity model data.
namespace ModelData
{
static double kappa_s = 1.0e6;
// for the block:
static double fixed_L = 0.25; // unit mm
// static double end_modulus_ratio = 0.4;
// static double end_coordinate = 30.0;
static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
static double slope = -2.6e-3; // slope for radius change 
static double r0 = 7.0e-2; // radius at root
//static double total_L = 30.0; // unit mm
// center: x=0.2; y=0.0

// Tether (penalty) force function for the solid block.
void block_tether_force_function(VectorValue<double>& F,
                                 const TensorValue<double>& /*FF*/,
                                 const libMesh::Point& X,
                                 const libMesh::Point& s,
                                 Elem* const /*elem*/,
                                 const vector<NumericVector<double>*>& /*system_data*/,
                                 double /*time*/,
                                 void* /*ctx*/)
{
    std::cout << "called block foce, check the code" << std::endl;
    F = kappa_s * (s - X);
    return;
} // block_tether_force_function

// Tether (penalty) force function for the thin beam.
void beam_tether_force_function(VectorValue<double>& F,
                                const TensorValue<double>& /*FF*/,
                                const libMesh::Point& X,
                                const libMesh::Point& s,
                                Elem* const /*elem*/,
                                const vector<NumericVector<double>*>& /*system_data*/,
                                double /*time*/,
                                void* /*ctx*/)
{
    const double r = sqrt((s(0) - 0) * (s(0) - 0));
    if (r <= fixed_L)
    {
        F = kappa_s * (s - X);
    }
    else
    {
        F.zero();
    }
    return;
} // beam_tether_force_function
// RE compute the arc-length to determine the slope
static bool Is_curved = false;
static double  delta_h = 0.005;
static double curved_a2 = 0.0; // second-order coefficient
static double curved_a1 = 0.0; // first-order coefficient
static double Beam_Length = 100.0;
double get_arc_length(const double x0)
{
  if (!Is_curved)
  return x0;
  // else
  int num_h = ceil(x0/delta_h);
  double sum_length =0.0;
  for (int k=0; k<num_h; k++)
  { double cur_x = k*delta_h;
    double dl_x = sqrt(1.0 + (2 * curved_a2 * cur_x + curved_a1) * (2 * curved_a2 * cur_x + curved_a1));
    sum_length = sum_length + dl_x * delta_h;
  }
  if (sum_length - 1.01* Beam_Length >= 0.0)
  pout << "check get_arc_length:" << sum_length << " > 1.01 Beam_L:" << Beam_Length << endl;
  return sum_length;

   
}

//
// Stress tensor function for the thin beam.
static double mu_s, lambda_s;
void beam_PK1_stress_function(TensorValue<double>& PP,
                              const TensorValue<double>& FF,
                              const libMesh::Point& /*X*/,
                              const libMesh::Point& s,
                              Elem* const /*elem*/,
                              const vector<NumericVector<double>*>& /*system_data*/,
                              double /*time*/,
                              void* /*ctx*/)
{      const double dist2root = sqrt((s(0) - 0) * (s(0) - 0));
    if (dist2root <= fixed_L) // we use penalty method to fix the root
    {
        PP = 0.0 * FF; // so we no longer compute stress here
	return ;
    }

 
    // 1)  compute the radius r(x) : r_x
    double arc_s = get_arc_length(s(0));
    const double r_x = slope * arc_s + r0;
    // 2) compute ratio_moduli;
    const double ratio_radius = r_x /r0; 
    const double ratio_moduli = ratio_radius *  ratio_radius*  ratio_radius*  ratio_radius;

   
        
        const TensorValue<double> CC = FF.transpose() * FF;
        const TensorValue<double> EE = 0.5 * (CC - II);
        const TensorValue<double> SS = lambda_s * EE.tr() * II + 2.0 * mu_s * EE;
        PP = ratio_moduli * FF * SS;
 
    return;
} // beam_PK1_stress_function


}
using namespace ModelData;

libMesh::Point tip_center;
// Function prototypes
static ofstream drag_stream, lift_stream, A_x_posn_stream, A_y_posn_stream, moment_stream;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& beam_mesh,
                      EquationSystems* beam_equation_systems,
                     //  Mesh& block_mesh,
                     //  EquationSystems* block_equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();

        const string beam_exodus_filename = app_initializer->getExodusIIFilename("beam");

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

               // >> obtain restart directory and restart number, added by walter
        const string restart_directory =app_initializer->getThisRestartDirectory();
        const int restart_number = app_initializer->getThisRestartNumber();

        // Create a simple FE mesh.
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;

     
	const double beam_h = input_db->getDouble("Beam_H");
	const double beam_l = input_db->getDouble("Beam_L");
	Is_curved = input_db->getBoolWithDefault("USING_INTRINSIC_CURVATURE", false);
	delta_h = ds * 0.05;
	curved_a2 = input_db->getDoubleWithDefault("CURVATURE_Coefficient_SecondOrder", 0.0);
 // second-order coefficient
	curved_a1 = input_db->getDoubleWithDefault("CURVATURE_Coefficient_FirstOrder", 0.0); // first-order coefficient
	Beam_Length = beam_l;
  if (Is_curved)
  { pout<<"+++++++++++++[Arc length mapping for stress]: delta_h=" <<delta_h << "; solid_mesh-size ds = " << ds << endl;
    pout<< "+++++++++++++[Arc length mapping for stress]: ( curved_a2, curved_a1) = " << "( "   << curved_a2 <<", " << curved_a1 << " ); beam_length = " << Beam_Length << endl;
  }
  else
  pout<<"++++++++++++++[Arch length mapping is not used]: no curvature" << endl;

        Mesh beam_mesh(NDIM);

// >> add 3D mesh (cylinder with radius) -- 02/03/2016 by walter

	if (NDIM==2)
{
		pout << "++++++++++++++[Geometry]:  we use 2D mesh ++++++++++++++++++++++ \n";
        string beam_elem_type = input_db->getString("BEAM_ELEM_TYPE");
		const double R = beam_h/2.0;
		//const double ds =input_db->getDouble("AXIAL_MESH_SIZE");
		
	    const int num_axial_elements = ceil(beam_l / ds); // elements
        MeshTools::Generation::build_square(beam_mesh,
                                            ceil(beam_l / ds),
                                            ceil(beam_h / ds),
                                            0.0,
                                            beam_l, //0.6,
                                            -0.5 * beam_h,
                                            0.5 * beam_h,
                                            Utility::string_to_enum<ElemType>(beam_elem_type));
	 const double R_tip = input_db->getDoubleWithDefault("TIP_RADIUS", R);
		if (R_tip == R)
			pout << "++++++++++++++[3D MESH]:  No taper, as R_tip = R ="<< R_tip <<" ++++++++++++++++++++++ \n";
		else if (R_tip < R)
			pout << "++++++++++++++[3D MESH]:  With Taper, as R_tip ="<< R_tip <<" R =" << R <<" ++++++++++++++++++++++ \n";
		else 
			pout << "++++++++++++++[3D MESH]:  Error: Wrong Taper: R_tip > R, as R_tip ="<< R_tip <<" R =" << R <<" ++++++++++++++++++++++ \n";
	
		double tip_scale = (R_tip / R);
	
	// deal with intrinsic curvature
		const bool is_with_curvature= input_db->getBoolWithDefault("USING_INTRINSIC_CURVATURE", false);
		const double coef_a2 = input_db->getDoubleWithDefault("CURVATURE_Coefficient_SecondOrder", 0.0);
		const double coef_a1 = input_db->getDoubleWithDefault("CURVATURE_Coefficient_FirstOrder", 0.0);
		int num_segments = input_db->getIntegerWithDefault("CURVATURE_Mapping_Number_of_segments",1); // number of segments when we compute the mapping
		if (is_with_curvature)
		{
		pout << "++++++++++++++[2D MESH]:  Intrinsic Curvature is used, (x1,0) ->(x2, f(x2)) ++++++++++++++++++++++ \n";
		pout << "++++++++++++++[2D MESH]:  f(x) = a2 x^2 + a1 x, where, a2 = " << coef_a2 << ", a1 = "<< coef_a1 << "++++++++++++++++++++++ \n";
		pout << "++++++++++++++[2D MESH]:  We do the mapping with number of segments as" << num_segments << "++++++++++++++++++++++ \n";
		}
		else
		{
		pout << "++++++++++++++[2D MESH]:  Intrinsic Curvature is not used ++++++++++++++++++++++ \n";
		num_segments = 1;
		pout << "++++++++++++++[2D MESH]:  We set number of segments as " << num_segments << "++++++++++++++++++++++ \n";
		}
		
	
		vector<double> vec_x2_at_nodes(num_axial_elements +1);
		const double exact_dz = beam_l / num_axial_elements;
		
		vec_x2_at_nodes[0] = 0.0; // (x1,0) ->(x2, f(x2)), when x1=0, we know x2 = 0

		tip_center=libMesh::Point(beam_l, 0.0, 0);
		double last_x1= 0.0;
		if(is_with_curvature)
		{  
			const double dseg = beam_l / num_segments;
			
			for (unsigned int kseg = 0; kseg < num_segments; kseg++)
				{   double cur_x2= dseg * (kseg +1); // x2 =  ds * (k+1)
					double g_x2 = sqrt(1.0 + (coef_a1 + 2.0 * coef_a2 * cur_x2) * (coef_a1 + 2.0 * coef_a2 * cur_x2)); //  g(x2) = \sqrt[1+ (a1 + 2 a2 x2)^2 ]
					double new_x1 = last_x1 + g_x2 * dseg; //x1[kseg+1] = x1[kseg] + g_x2 * dseg;
					unsigned int new_ind = ceil(new_x1 / exact_dz);
					unsigned int last_ind = ceil(last_x1 / exact_dz);
					if ( (last_ind < new_ind) && (last_ind < num_axial_elements +1) ) // this interval contains one node: vec_x2_at_nodes[last_ind]
					{   double last_x2 = kseg * dseg; 
						double new_x2 = (kseg+1) * dseg;
						double cur_x1 = last_ind * exact_dz;
						vec_x2_at_nodes[last_ind] = (cur_x1 - last_x1) / (new_x1 - last_x1) * new_x2 + (new_x1 - cur_x1) / (new_x1 - last_x1) * last_x2;
						
					}
					// next step
					last_x1 = new_x1;
					 
				}
		// we check the vec_x2(x1);
			
		pout << "++++++++++++++[2D MESH]:  finish the x2-x1 mapping with segments as " << num_segments << "++++++++++++++++++++++ \n";	
		// get the new tip_center
				double cur_x2 = vec_x2_at_nodes[num_axial_elements];				
				double dfdx= coef_a1 + 2.0* coef_a2 * cur_x2;
				double fx= coef_a1 * cur_x2 + coef_a2 * cur_x2 * cur_x2;
				double g_x2 = sqrt(1.0 + dfdx * dfdx); // g(x2)

				tip_center(0) = cur_x2;
				tip_center(1) = fx;
				tip_center(2) = 0;	
		pout << "++++++++++++++[2D MESH]:  Now the coordinates of tip-center are: " << tip_center<< "++++++++++++++++++++++ \n";	
		
		}
		// change the node_coordinates
		Mesh::node_iterator       it_nd      = beam_mesh.nodes_begin(); //mesh.active_local_elements_begin();//mesh.elements_begin();
    	        const Mesh::node_iterator it_last_nd = beam_mesh.nodes_end(); //mesh.active_local_elements_end();//mesh.elements_end();
    	for ( ; it_nd != it_last_nd ; ++it_nd) {
        	Node* node = *it_nd;
		// step 1) new_x = old_z; new_y = old_x; new_z = old_y
			double new_x = (*node)(0);
			double new_y = (*node)(1);
			

		// step 2) add varing radius for taper
			double current_scale = (new_x / beam_l) * tip_scale + ( 1.0 - new_x / beam_l) * 1.0;
			new_y = new_y* current_scale;
				
		// step 3) add intrinsic curvature (see the note)
			if (is_with_curvature)
			{
				 // pout << "++++++++++++++[3D MESH]:  Begin the curvature mapping ++++++++++++++++++++++ \n";		
				
				// (x1,y1,z1) --mapped to-->(x2,y2,z1)
				// first, we need to compute x2
				unsigned int cur_ind = ceil ( new_x / exact_dz + 0.2 )  - 1; // to get the index;
				double cur_x2 = vec_x2_at_nodes[cur_ind];				
				double dfdx= coef_a1 + 2.0* coef_a2 * cur_x2;
				double fx= coef_a1 * cur_x2 + coef_a2 * cur_x2 * cur_x2;
				double g_x2 = sqrt(1.0 + dfdx * dfdx); // g(x2)

				(*node)(0) = cur_x2 - dfdx * new_y / g_x2;
				(*node)(1) = fx + new_y /g_x2 ;
								
				/*				
				pout <<" cur_ind is = " << cur_x2 <<"; new_x2 = " << cur_x2 - dfdx * new_y / g_x2 << endl;
				pout <<"previous node cords: (" << new_x <<", " << new_y <<", " << new_z << endl;
				pout <<"current node cords: (" << (*node)(0) << ", " << (*node)(1)<<", " << (*node)(2)<< endl;
				*/
								
			}
			else
			{
				(*node)(0) = new_x;
				(*node)(1) = new_y;
				
			}
		
		
	   } // for 
        
}
	else // NDIM ==3
{
		pout << "++++++++++++++[Geometry]:  we use 3D mesh ++++++++++++++++++++++ \n";
        Mesh block_mesh(2); // circle mesh
		string block_elem_type = input_db->getString("CSA_ELEM_TYPE"); // currently, tri3 or tri6
		const double dr = input_db->getDouble("CSA_MESH_SIZE");
		const double R = input_db->getDouble("ROOT_RADIUS");
		const double dz =input_db->getDouble("AXIAL_MESH_SIZE");
		const int num_circum_nodes = ceil(2.0 * M_PI * R / dr);
	    const int num_axial_elements = ceil(beam_l / dz); // elements
		pout << "++++++++++++++[Geometry]:  First is 2D CSA mesh ++++++++++++++++++++++ \n";
		pout << "+++++++++++++++++++++++++++ [2D CSA]: (R, dz) = (" << R << ", " <<dz << ")++++++++++++++++++++++ \n";
		
        for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                block_mesh.add_point(libMesh::Point( R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(block_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(block_elem_type);
            triangle.desired_area() = sqrt(3.0) / 4.0 * dr * dr;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
            block_mesh.prepare_for_use();
		pout << "++++++++++++++[Geometry]:  finish generatign 2D mesh ++++++++++++++++++++++ \n";
		const RealVectorValue extrusion_vector (0.0, 0.0, beam_l);
		pout << "++++++++++++++[Geometry]:  begin to deal with 3D mesh ++++++++++++++++++++++ \n";
		const double R_tip = input_db->getDoubleWithDefault("TIP_RADIUS", R);
		if (R_tip == R)
			pout << "++++++++++++++[3D MESH]:  No taper, as R_tip = R ="<< R_tip <<" ++++++++++++++++++++++ \n";
		else if (R_tip < R)
			pout << "++++++++++++++[3D MESH]:  With Taper, as R_tip ="<< R_tip <<" R =" << R <<" ++++++++++++++++++++++ \n";
		else 
			pout << "++++++++++++++[3D MESH]:  Error: Wrong Taper: R_tip > R, as R_tip ="<< R_tip <<" R =" << R <<" ++++++++++++++++++++++ \n";
	
		double tip_scale = (R_tip / R);
		
		// deal with intrinsic curvature
		const bool is_with_curvature= input_db->getBoolWithDefault("USING_INTRINSIC_CURVATURE", false);
		const double coef_a2 = input_db->getDoubleWithDefault("CURVATURE_Coefficient_SecondOrder", 0.0);
		const double coef_a1 = input_db->getDoubleWithDefault("CURVATURE_Coefficient_FirstOrder", 0.0);
		int num_segments = input_db->getIntegerWithDefault("CURVATURE_Mapping_Number_of_segments",1); // number of segments when we compute the mapping
		if (is_with_curvature)
		{
		pout << "++++++++++++++[3D MESH]:  Intrinsic Curvature is used, (x1,0) ->(x2, f(x2)) ++++++++++++++++++++++ \n";
		pout << "++++++++++++++[3D MESH]:  f(x) = a2 x^2 + a1 x, where, a2 = " << coef_a2 << ", a1 = "<< coef_a1 << "++++++++++++++++++++++ \n";
		pout << "++++++++++++++[3D MESH]:  We do the mapping with number of segments as" << num_segments << "++++++++++++++++++++++ \n";
		}
		else
		{
		pout << "++++++++++++++[3D MESH]:  Intrinsic Curvature is not used ++++++++++++++++++++++ \n";
		num_segments = 1;
		pout << "++++++++++++++[3D MESH]:  We set number of segments as " << num_segments << "++++++++++++++++++++++ \n";
		}
		
	
		vector<double> vec_x2_at_nodes(num_axial_elements +1);
		const double exact_dz = beam_l / num_axial_elements;
		
		vec_x2_at_nodes[0] = 0.0; // (x1,0) ->(x2, f(x2)), when x1=0, we know x2 = 0

		tip_center=libMesh::Point(beam_l, 0.0, 0);
		double last_x1= 0.0;
		if(is_with_curvature)
		{  
			const double dseg = beam_l / num_segments;
			
			for (unsigned int kseg = 0; kseg < num_segments; kseg++)
				{   double cur_x2= dseg * (kseg +1); // x2 =  ds * (k+1)
					double g_x2 = sqrt(1.0 + (coef_a1 + 2.0 * coef_a2 * cur_x2) * (coef_a1 + 2.0 * coef_a2 * cur_x2)); //  g(x2) = \sqrt[1+ (a1 + 2 a2 x2)^2 ]
					double new_x1 = last_x1 + g_x2 * dseg; //x1[kseg+1] = x1[kseg] + g_x2 * dseg;
					unsigned int new_ind = ceil(new_x1 / exact_dz);
					unsigned int last_ind = ceil(last_x1 / exact_dz);
					if ( (last_ind < new_ind) && (last_ind < num_axial_elements +1) ) // this interval contains one node: vec_x2_at_nodes[last_ind]
					{   double last_x2 = kseg * dseg; 
						double new_x2 = (kseg+1) * dseg;
						double cur_x1 = last_ind * exact_dz;
						vec_x2_at_nodes[last_ind] = (cur_x1 - last_x1) / (new_x1 - last_x1) * new_x2 + (new_x1 - cur_x1) / (new_x1 - last_x1) * last_x2;
						
					}
					// next step
					last_x1 = new_x1;
					 
				}
		// we check the vec_x2(x1);
			
		pout << "++++++++++++++[3D MESH]:  finish the x2-x1 mapping with segments as " << num_segments << "++++++++++++++++++++++ \n";	
		// get the new tip_center
				double cur_x2 = vec_x2_at_nodes[num_axial_elements];				
				double dfdx= coef_a1 + 2.0* coef_a2 * cur_x2;
				double fx= coef_a1 * cur_x2 + coef_a2 * cur_x2 * cur_x2;
				double g_x2 = sqrt(1.0 + dfdx * dfdx); // g(x2)

				tip_center(0) = cur_x2;
				tip_center(1) = fx;
				tip_center(2) = 0;	
		pout << "++++++++++++++[3D MESH]:  Now the coordinates of tip-center are: " << tip_center<< "++++++++++++++++++++++ \n";	
		
		}

        MeshTools::Generation::build_extrusion(beam_mesh,
                                                block_mesh,num_axial_elements, extrusion_vector);
		// change the node_coordinates
		Mesh::node_iterator       it_nd      = beam_mesh.nodes_begin(); //mesh.active_local_elements_begin();//mesh.elements_begin();
    	const Mesh::node_iterator it_last_nd = beam_mesh.nodes_end(); //mesh.active_local_elements_end();//mesh.elements_end();
    	for ( ; it_nd != it_last_nd ; ++it_nd) {
        	Node* node = *it_nd;
		// step 1) new_x = old_z; new_y = old_x; new_z = old_y
			double new_x = (*node)(2);
			double new_y = (*node)(0);
			double new_z = (*node)(1);

		// step 2) add varing radius for taper
			double current_scale = (new_x / beam_l) * tip_scale + ( 1.0 - new_x / beam_l) * 1.0;
			new_y = new_y* current_scale;
			new_z = new_z* current_scale;	
		// step 3) add intrinsic curvature (see the note)
			if (is_with_curvature)
			{
				 // pout << "++++++++++++++[3D MESH]:  Begin the curvature mapping ++++++++++++++++++++++ \n";		
				
				// (x1,y1,z1) --mapped to-->(x2,y2,z1)
				// first, we need to compute x2
				unsigned int cur_ind = ceil ( new_x / exact_dz + 0.2 )  - 1; // to get the index;
				double cur_x2 = vec_x2_at_nodes[cur_ind];				
				double dfdx= coef_a1 + 2.0* coef_a2 * cur_x2;
				double fx= coef_a1 * cur_x2 + coef_a2 * cur_x2 * cur_x2;
				double g_x2 = sqrt(1.0 + dfdx * dfdx); // g(x2)

				(*node)(0) = cur_x2 - dfdx * new_y / g_x2;
				(*node)(1) = fx + new_y /g_x2 ;
				(*node)(2) = new_z;				
				/*				
				pout <<" cur_ind is = " << cur_x2 <<"; new_x2 = " << cur_x2 - dfdx * new_y / g_x2 << endl;
				pout <<"previous node cords: (" << new_x <<", " << new_y <<", " << new_z << endl;
				pout <<"current node cords: (" << (*node)(0) << ", " << (*node)(1)<<", " << (*node)(2)<< endl;
				*/
								
			}
			else
			{
				(*node)(0) = new_x;
				(*node)(1) = new_y;
				(*node)(2) = new_z;
			}
		
		
	   } // for 


}

// << add 3D mesh (cylinder with radius) -- 02/03/2016 by walter
		beam_mesh.prepare_for_use();
		pout << "++++++++++++++[Geometry]: finish generating the beam mesh ++++++++++++++++++++++ \n";
        
        fixed_L = input_db->getDouble("FIXED_L");
        mu_s = input_db->getDouble("MU_S");
        lambda_s = input_db->getDouble("LAMBDA_S");
        kappa_s = input_db->getDouble("KAPPA_S");
		slope =input_db->getDouble("SLOPE_RADIUS");
		r0 = input_db->getDouble("ROOT_RADIUS");
	// >> begin to print out
        pout << "\n";
		pout << "Modulus = " << mu_s <<"; Bulk Modulus =" << lambda_s << "\n";
		pout << "SLOPE of radius change: " << slope << "; Radius at root " << r0 << "\n";
	// << end print out
	// end_modulus_ratio = input_db->getDouble("END_MODULUS_RATIO"); // not used, change on 01/29/2016
	// end_coordinate = input_db->getDouble("END_COORDINATE"); // not used change on 01/29/2016
        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &beam_mesh, 
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels")
			   , restart_directory,  restart_number
			   );
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);


        // Configure the IBFE solver.
	

        IBFEMethod::LagBodyForceFcnData beam_tether_force_data(beam_tether_force_function);
        IBFEMethod::PK1StressFcnData beam_PK1_stress_data(beam_PK1_stress_function);
        ib_method_ops->registerLagBodyForceFunction(beam_tether_force_data);
        ib_method_ops->registerPK1StressFunction(beam_PK1_stress_data);
 
        EquationSystems* beam_equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();
	
        FEDataManager* fe_data_manager = ib_method_ops->getFEDataManager(0);
        Pointer<IBFEPostProcessor> ib_post_processor = new IBFECentroidPostProcessor("IBFEPostProcessor", fe_data_manager);

        ib_post_processor->registerTensorVariable(
            "FF", MONOMIAL, CONSTANT, IBFEPostProcessor::FF_fcn);

        std::pair<IBTK::TensorMeshFcnPtr,void*> PK1_dev_stress_fcn_data(beam_PK1_stress_function,static_cast<void*>(NULL));
        ib_post_processor->registerTensorVariable(
            "sigma_dev", MONOMIAL, CONSTANT,
            IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
            std::vector<unsigned int>(), &PK1_dev_stress_fcn_data);
	// >> add pressure interpolation
	Pointer<IBFEPostProcessor> ib_pressure_processor =
            new IBFECentroidPostProcessor("IBFEPressureProcessor", fe_data_manager);
	 std::vector<double> vec_ones (3,1.0);
        HierarchyGhostCellInterpolation::InterpolationTransactionComponent p_ghostfill(
/*data_idx*/ -1, "LINEAR_REFINE", /*use_cf_bdry_interpolation*/ false, "CONSERVATIVE_COARSEN", "LINEAR");
	    
		FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR", QGAUSS, FIFTH, /*use_adaptive_quadrature*/ false,
/*point_density*/ 2.0, /*use_consistent_mass_matrix*/ true, false, vec_ones);



	ib_pressure_processor->registerInterpolatedScalarEulerianVariable("pressure_f",
                                                                      MONOMIAL, //LAGRANGE,
                                                                      CONSTANT, // FIRST,
                                                                      navier_stokes_integrator->getPressureVariable(),
                                                                      navier_stokes_integrator->getCurrentContext(),
								      p_ghostfill, p_interp_spec);   
	// << add pressure interpolation


        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        // AutoPtr<ExodusII_IO> block_exodus_io(uses_exodus ? new ExodusII_IO(block_mesh) : NULL);
        AutoPtr<ExodusII_IO> beam_exodus_io(uses_exodus ? new ExodusII_IO(beam_mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
	ib_post_processor->initializeFEData();
	ib_pressure_processor->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
		ib_post_processor->postProcessData(loop_time); 
		ib_pressure_processor->postProcessData(loop_time);
                beam_exodus_io->write_timestep(
                    beam_exodus_filename, *beam_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            
	                   std::ostringstream file_name;
            file_name << beam_exodus_filename +"_"
                      << std::setw(6)
                      << std::setfill('0')
                      << std::right
                      << iteration_num;
	    GMVIO(beam_mesh).write_equation_systems(file_name.str()+".gmv",*beam_equation_systems);
            }
        }

        // Open streams to save lift and drag coefficients.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.open("C_D.curve", ios_base::out | ios_base::trunc);
            lift_stream.open("C_L.curve", ios_base::out | ios_base::trunc);
	    moment_stream.open("Moment.curve", ios_base::out | ios_base::trunc);
            A_x_posn_stream.open("A_x.curve", ios_base::out | ios_base::trunc);
            A_y_posn_stream.open("A_y.curve", ios_base::out | ios_base::trunc);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            time_integrator->advanceHierarchy(dt);
            loop_time += dt;
	    
	  //  ib_pressure_processor->postProcessData(loop_time);
            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
	    

	    
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
			   ib_pressure_processor->postProcessData(loop_time);
             		ib_post_processor->postProcessData(loop_time);     
		     		beam_exodus_io->write_timestep(
                        beam_exodus_filename, *beam_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                
                    std::ostringstream file_name;
            		file_name << beam_exodus_filename +"_"
                      << std::setw(6)
                      << std::setfill('0')
                      << std::right
                      << iteration_num;
            		GMVIO(beam_mesh).write_equation_systems(file_name.str()+".gmv",*beam_equation_systems);


				}
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
           	    ib_method_ops->writeRestartEquationSystems(restart_dump_dirname, iteration_num);
	    }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                pout << "\nWriting state data...\n\n";
                postprocess_data(patch_hierarchy,
                                 navier_stokes_integrator,
                                 beam_mesh,
                                 beam_equation_systems,
                                 //,  block_mesh,
                                 //, // block_equation_systems,
                                 iteration_num,
                                 loop_time,
                                 postproc_data_dump_dirname);
            }
        }

        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            drag_stream.close();
            lift_stream.close();
	    moment_stream.close();
            A_x_posn_stream.close();
            A_y_posn_stream.close();
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return 0;
} // main

void postprocess_data(Pointer<PatchHierarchy<NDIM> > /*patch_hierarchy*/,
                      Pointer<INSHierarchyIntegrator> /*navier_stokes_integrator*/,
                      Mesh& beam_mesh,
                      EquationSystems* beam_equation_systems,
                      const int /*iteration_num*/,
                      const double loop_time,
                      const string& /*data_dump_dirname*/)
{
    double F_integral[NDIM];
    // >> add vector for moment
    double M_integral=0.0;
    // <<
    for (unsigned int d = 0; d < NDIM; ++d) F_integral[d] = 0.0;
    Mesh* mesh[1] = { &beam_mesh };
    EquationSystems* equation_systems[1] = { beam_equation_systems};
    for (unsigned int k = 0; k < 1; ++k)
    {
        System& F_system = equation_systems[k]->get_system<System>(IBFEMethod::FORCE_SYSTEM_NAME);
        NumericVector<double>* F_vec = F_system.solution.get();
        NumericVector<double>* F_ghost_vec = F_system.current_local_solution.get();
        F_vec->localize(*F_ghost_vec);
        DofMap& F_dof_map = F_system.get_dof_map();
        std::vector<std::vector<unsigned int> > F_dof_indices(NDIM);
        AutoPtr<FEBase> fe(FEBase::build(NDIM, F_dof_map.variable_type(0)));
        AutoPtr<QBase> qrule = QBase::build(QGAUSS, NDIM, FIFTH);
        fe->attach_quadrature_rule(qrule.get());
        const std::vector<std::vector<double> >& phi = fe->get_phi();
        const std::vector<double>& JxW = fe->get_JxW();
	//>> add coordinates to compute moment
	const std::vector<libMesh::Point>& xyz_qps = fe->get_xyz();
	
	//<<
        boost::multi_array<double, 2> F_node;
        const MeshBase::const_element_iterator el_begin = mesh[k]->active_local_elements_begin();
        const MeshBase::const_element_iterator el_end = mesh[k]->active_local_elements_end();
        for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
        {
            Elem* const elem = *el_it;
            fe->reinit(elem);
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                F_dof_map.dof_indices(elem, F_dof_indices[d], d);
            }
            const int n_qp = qrule->n_points();
            const int n_basis = F_dof_indices[0].size();
            get_values_for_interpolation(F_node, *F_ghost_vec, F_dof_indices);
            for (int qp = 0; qp < n_qp; ++qp)
            {
                for (int k = 0; k < n_basis; ++k)
                {
                    for (int d = 0; d < NDIM; ++d)
                    {
                        F_integral[d] += F_node[k][d] * phi[k][qp] * JxW[qp];
			
                    }
		    M_integral += (F_node[k][1] * xyz_qps[qp](0) -  F_node[k][0] * xyz_qps[qp](1)) * phi[k][qp] * JxW[qp];
                }
            }
        }
    }
    SAMRAI_MPI::sumReduction(F_integral, NDIM);
    SAMRAI_MPI::sumReduction(M_integral);
    if (SAMRAI_MPI::getRank() == 0)
    {
        drag_stream.precision(12);
        drag_stream.setf(ios::fixed, ios::floatfield);
        drag_stream << loop_time << " " << -F_integral[0] << endl;
        lift_stream.precision(12);
        lift_stream.setf(ios::fixed, ios::floatfield);
        lift_stream << loop_time << " " << -F_integral[1] << endl;
        moment_stream.precision(12);
	moment_stream.setf(ios::fixed, ios::floatfield);
        moment_stream << loop_time << " " << -M_integral << endl;
    }

    System& X_system = beam_equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    NumericVector<double>* X_vec = X_system.solution.get();
    AutoPtr<NumericVector<Number> > X_serial_vec = NumericVector<Number>::build(X_vec->comm());
    X_serial_vec->init(X_vec->size(), true, SERIAL);
    X_vec->localize(*X_serial_vec);
    DofMap& X_dof_map = X_system.get_dof_map();
    vector<unsigned int> vars(2);
    vars[0] = 0;
    vars[1] = 1;
    MeshFunction X_fcn(*beam_equation_systems, *X_serial_vec, X_dof_map, vars);
    X_fcn.init();
    DenseVector<double> X_A(2);
    X_fcn(tip_center, 0.0, X_A);
    if (SAMRAI_MPI::getRank() == 0)
    {
        A_x_posn_stream.precision(12);
        A_x_posn_stream.setf(ios::fixed, ios::floatfield);
        A_x_posn_stream << loop_time << " " << X_A(0) << endl;
        A_y_posn_stream.precision(12);
        A_y_posn_stream.setf(ios::fixed, ios::floatfield);
        A_y_posn_stream << loop_time << " " << X_A(1) << endl;
    }
    return;
} // postprocess_data
