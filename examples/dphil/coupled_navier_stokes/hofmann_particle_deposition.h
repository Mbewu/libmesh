
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
		HofmannParticleDeposition (EquationSystems& es_in, std::vector<std::vector<double> >& element_data_in, unsigned int _num_generations);

		// efficiency of deposition of particle being deposited at this time
		void calculate_deposition_efficiency(unsigned int t_step_in);
		// calculate the flow rate at all times for all branches and put into HofmannParticleDeposition::flow_function
		void calculate_flow_rate();

		void calculate_total_efficiency_function();

		std::vector<double> get_total_efficiency() { return total_efficiency_function;};

		std::vector<std::vector<double> > get_acinar_efficiency_per_generation() { return alveolar_efficiency_per_generation;};
		std::vector<std::vector<double> > get_tb_efficiency_per_generation() { return tb_efficiency_per_generation;};
		std::vector<double> get_non_deposited_particle_weight() { return non_deposited_particle_weight;};
		std::vector<double> get_output_particle_weight() { return output_particle_weight;};
		int get_total_particles_inhaled() { return total_particles_inhaled;}

	
	private:
		
		double calculate_sedimentation_efficiency(double pipe_angle_from_gravity,double pipe_radius,double pipe_length,double mean_pipe_velocity,double time_in_branch);
		double calculate_brownian_efficiency(double pipe_radius, double pipe_length,double mean_pipe_velocity,double time_in_branch);
		double calculate_impaction_efficiency(double branching_angle,double stokes_number, double reynolds_number);

		void construct_temporary_alveolated_1d_tree (unsigned int parent_id);

		double time_scaling(double time_in);

		// efficiency of deposition of particle being deposited at this time
		void calculate_deposition_efficiency_in_tb_daughters
		(unsigned int t_step_in,  unsigned int parent_id,double parent_end_time, double parent_end_flow, double parent_end_weight);

		//calculate the deposition efficiency in branch_id
		void calculate_deposition_efficiency_in_tb_parent
		(unsigned int t_step_in,  unsigned int branch_id,double daughter_end_time, double daughter_end_flow, double daughter_end_weight);

		//calculate deposition in daughters of parent_alveolar_id
		void calculate_deposition_efficiency_in_alveolar_daughters
		(unsigned int t_step_in,  unsigned int parent_alveolar_id,double parent_end_time, double parent_end_flow, double parent_end_weight);

		//calcualte deposition in branch alveolar_id
		void calculate_deposition_efficiency_in_alveolar_parent
		(unsigned int t_step_in,  unsigned int alveolar_id,double daughter_end_time, double daughter_end_flow, double daughter_end_weight);

		// calculate the flow rate at daughters of parent_id with parent_flow
		void calculate_flow_at_daughters(unsigned int t_step, unsigned int parent_id, double parent_flow);

		// calculate the efficiency on branch_id with flow rate branch_flow
		double calculate_efficiency(double branch_flow,double pipe_radius, double pipe_length, double gravity_angle,double branching_angle,double time_in_branch);

		// calculate the efficiency of deposition in the alveoli
		double calculate_deposition_in_alveoli(double time_spent_in_alveoli);

		// calculate the sedimentation efficiency in a sphere
		double calculate_sedimentation_efficiency_in_sphere(double time_spent_in_alveoli);
	
		// calculate the brownian efficiency in a sphere
		double calculate_brownian_efficiency_in_sphere(double time_spent_in_alveoli);


		//returns false if at end of time
		bool calculate_flow_at_time(double& flow,double time_in, unsigned int branch_no);




		EquationSystems* es;
		std::vector<std::vector<double> > element_data;	
		std::vector<std::vector<double> > efficiency_function; // [timestep][element_no], efficiency of deposition of particle being deposited
		std::vector<double > total_efficiency_function; // [timestep][element_no], efficiency of deposition of particle being deposited
		std::vector<std::vector<double> > alveolar_efficiency_per_generation;	// [timestep][generation] alveolar efficiency of deposition summed over all branches in this generation
		std::vector<std::vector<double> > tb_efficiency_per_generation;	// [timestep][generation] tracheo-bronchial efficiency of deposition summed over all branches in this generation
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
		std::vector<double> output_particle_weight;				// this is the particle weight that is exhaled
		std::vector<double> non_deposited_particle_weight;	// [timestep] this is the particle weight that wasn't deposited at all
		double alveolus_radius;
		double sedimentation_velocity;
		unsigned int num_generations;
		double dt;

		std::vector<std::vector<double> > temporary_alveolar_tree;	

};
#endif //__hofmann_particle_deposition_h__

