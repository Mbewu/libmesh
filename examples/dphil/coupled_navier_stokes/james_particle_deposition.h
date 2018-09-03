

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

#include "airway.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifndef __james_particle_deposition_h__
#define __james_particle_deposition_h__

// headers omitted for brevity
class JamesParticleDeposition
{
	public:
		JamesParticleDeposition ();

		// setup everything
		// edit: removed default argument for surface_boundaries that was an empty vector
		void setup(EquationSystems& es_in,EquationSystems& es_1d_in, std::vector<Airway>& _airway_data, std::vector<unsigned int> _subdomains_1d, std::vector<SurfaceBoundary* > surface_boundaries, EquationSystems& es_1d_terminal_in, std::vector<std::vector<Airway> >& _terminal_airway_data, std::vector<unsigned int> _terminal_airway_to_airway_map);

		// calculate the deposition probability of each airway
		void calculate_airway_deposition_probability();

		/* do later
		// calculate the deposition probability of each terminal airway
		void calculate_terminal_airway_deposition_probability();
		*/

		// based on the fraction exited 3d, or just start at zero
		void calculate_total_deposition_fraction(std::vector<double> total_fraction_exited_3d_surface=std::vector<double>(), std::vector<std::vector<unsigned int> > airway_elem_id_starts=std::vector<std::vector<unsigned int> >());

		// based on the fraction exited 3d, or just start at zero
		void calculate_time_step_deposition_fraction(std::vector<std::vector<unsigned int> > airway_elem_id_starts, double _time, double _dt);

		// to initialise particle fraction before first time step, assumes 1 inlet and t=0
		void add_particle_fractions_to_tree(std::vector<double> total_fraction_exited_3d_surface, std::vector<std::vector<unsigned int> > airway_elem_id_starts, double _time);

		// return the total particle fraction in the airways
		double get_total_particle_fraction() { return total_particle_fraction; };

		// get the total deposition fraction and exit fraction
		double get_total_deposition_fraction() { return total_deposition_fraction; };
		double get_total_deposition_fraction_sed() { return total_deposition_fraction_sed; };
		double get_total_deposition_fraction_imp() { return total_deposition_fraction_imp; };
		double get_total_deposition_fraction_dif() { return total_deposition_fraction_dif; };
		double get_total_exit_fraction() { return total_exit_fraction; };
		double get_total_deposition_fraction_this_time_step() { return total_deposition_fraction_this_time_step; };
		double get_total_exit_fraction_this_time_step() { return total_exit_fraction_this_time_step; };


	
	private:
		
		double calculate_sedimentation_probability(double pipe_angle_from_gravity,double pipe_radius,double pipe_length,double mean_pipe_velocity,double time_in_branch);
		double calculate_brownian_probability(double pipe_radius, double pipe_length,double mean_pipe_velocity,double time_in_branch);
		double calculate_impaction_probability(double branching_angle,double stokes_number, double reynolds_number);


		void calculate_airway_total_deposition_fraction(unsigned int element_number, double entering_particle_density);

		void calculate_airway_time_step_deposition_fraction(unsigned int element_number);

		// return the airway deposition prob
		void get_airway_deposition_probability(unsigned int element_number, double mean_velocity, double& p_airway_total, double& p_airway_sed, double& p_airway_imp, double& p_airway_dif);

		EquationSystems* es;
		EquationSystems* es_1d;
		EquationSystems* es_1d_terminal;
		std::vector<Airway>* p_airway_data;	
		std::vector<std::vector<Airway> >* p_terminal_airway_data;	
		double particle_density;
		double cunningham_correction_factor;
		double particle_radius;
		double particle_diameter;
		double gravity_magnitude;
		double density;
		double viscosity;
		double diffusion_coefficient;
		double sedimentation_velocity;

		double time;
		double dt;

		std::vector<double> branching_angles;	// [element_no]
		std::vector<double> gravity_angles;	// [element_no]
		std::vector<double> airway_deposition_probability;	// [element_no]
		std::vector<double> airway_deposition_probability_sed;	// [element_no]
		std::vector<double> airway_deposition_probability_imp;	// [element_no]
		std::vector<double> airway_deposition_probability_dif;	// [element_no]
		std::vector<double> airway_flow_rate;	// [element_no]

		double total_deposition_fraction_this_time_step;
		double total_exit_fraction_this_time_step;
		double total_deposition_fraction;
		double total_deposition_fraction_sed;
		double total_deposition_fraction_imp;
		double total_deposition_fraction_dif;
		double total_exit_fraction;
		double total_particle_fraction;

		std::vector<unsigned int> subdomains_1d;

		bool deposition_verbose;


		// terminal data
		std::vector<double> terminal_branching_angles;	// [element_no]
		std::vector<double> terminal_gravity_angles;	// [element_no]
		bool alveolar_deposition;
		std::vector<std::vector<double> > terminal_airway_deposition_probability;	// [element_no]
		std::vector<std::vector<double> > terminal_airway_deposition_probability_sed;	// [element_no]
		std::vector<std::vector<double> > terminal_airway_deposition_probability_imp;	// [element_no]
		std::vector<std::vector<double> > terminal_airway_deposition_probability_dif;	// [element_no]
		std::vector<std::vector<double> > terminal_airway_flow_rate;	// [element_no]

		double total_terminal_deposition_fraction_this_time_step;
		double total_terminal_exit_fraction_this_time_step;
		double total_terminal_deposition_fraction;
		double total_terminal_deposition_fraction_sed;
		double total_terminal_deposition_fraction_imp;
		double total_terminal_deposition_fraction_dif;
		double total_terminal_exit_fraction;
		double total_terminal_particle_fraction;

		std::vector<unsigned int> terminal_airway_to_airway_map;



};
#endif //__hofmann_particle_deposition_h__

