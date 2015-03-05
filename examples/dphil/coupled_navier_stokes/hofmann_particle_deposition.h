
// C++ include files that we may need
#include <iostream>
#include <algorithm>
#include <sstream>
#include <math.h>


//may need this
#include "libmesh/libmesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/perf_log.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h"

// For systems of equations the \p DenseSubMatrix
// and \p DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

// The definition of a geometric element
#include "libmesh/elem.h"

//for parsing command line nicely
#include "libmesh/getpot.h"

//include for general assembly class
#include "libmesh/system.h"

#include "surface_boundary.h"


#include "libmesh/point_locator_base.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifndef __hofmann_particle_deposition_h__
#define __hofmann_particle_deposition_h__

// headers omitted for brevity
class HofmannParticleDeposition
{
	public:
		HofmannParticleDeposition (EquationSystems& es_in, std::vector<std::vector<double> >& element_data_in);

		// efficiency of deposition of particle being deposited at this time
		void calculate_deposition_efficiency(unsigned int t_step_in);
		// efficiency of deposition of particle being deposited at this time
		void calculate_deposition_efficiency_at_daughters(unsigned int t_step_in,  unsigned int parent_id,double parent_end_time, double parent_end_weight, double parent_end_particle_weight);
		// calculate the flow rate at all times for all branches and put into HofmannParticleDeposition::flow_function
		void calculate_flow_rate();
		// calculate the flow rate at daughters of parent_id with parent_flow
		void calculate_flow_at_daughters(unsigned int t_step, unsigned int parent_id, double parent_flow);
		// calculate the efficiency on branch_id with flow rate branch_flow
		double calculate_efficiency(unsigned int branch_id,double branch_flow);

		void calculate_total_efficiency_function();

		std::vector<double> get_total_efficiency() { return total_efficiency_function;};

	
	private:
		
		double calculate_sedimentation_efficiency(double pipe_angle_from_gravity,double pipe_radius,double pipe_length,double mean_pipe_velocity);
		double calculate_brownian_efficiency(double pipe_radius, double pipe_length,double mean_pipe_velocity);
		double calculate_impaction_efficiency(double branching_angle,double stokes_number);

		double time_scaling(double time_in);
		//returns false if at end of time
		bool calculate_flow_at_time(double& flow,double time_in, unsigned int branch_no);

		EquationSystems* es;
		std::vector<std::vector<double> > element_data;	
		std::vector<std::vector<double> > efficiency_function; // [timestep][element_no], efficiency of deposition of particle being deposited
		std::vector<double > total_efficiency_function; // [timestep][element_no], efficiency of deposition of particle being deposited
		std::vector<std::vector<double> > flow_function;	// [timestep][element_no]
		std::vector<double> branching_angles;	// [element_no]
		std::vector<double> gravity_angles;	// [element_no]
		std::vector<double> times;	// [timestep][element_no]
		double particle_density;
		double cunningham_correction_factor;
		double particle_radius;
		double gravity_magnitude;
		double viscosity;
		double diffusion_coefficient;
		int total_particles_inhaled;
		bool end_sim_at_end_time;

};
#endif //__hofmann_particle_deposition_h__

