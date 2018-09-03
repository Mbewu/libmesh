
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

#ifndef __airway_h__
#define __airway_h__

// headers omitted for brevity
class Airway
{
	public:
		Airway ();

		void set_generation(int _generation) { generation = _generation; };
		int get_generation();

		void set_is_daughter_1(bool _is_daughter_1) { is_daughter_1 = _is_daughter_1; };
		bool get_is_daughter_1() { return is_daughter_1;};

		void add_sibling(unsigned int sibling) { siblings.push_back(sibling); };
		std::vector<unsigned int> get_siblings() { return siblings; };

		unsigned int get_num_siblings() { return siblings.size(); };

		void set_primary_sibling(int _primary_sibling) { primary_sibling = _primary_sibling; };
		int get_primary_sibling();
		bool has_primary_sibling();

		void set_length(double _length) { length = _length; };
		double get_length();

		void set_radius(double _radius) { radius = _radius; };
		double get_radius();

		double get_area() { return M_PI*radius*radius; };

		void set_order(int _order) { order = _order; };
		int get_order();

		void set_flow_rate(double _flow_rate) { flow_rate = _flow_rate; };
		double get_flow_rate() { return flow_rate; };

		void set_efficiency(double _efficiency) { efficiency = _efficiency; };
		void add_efficiency(double _efficiency) { efficiency += _efficiency; };
		double get_efficiency() { return efficiency; };

		void set_num_alveolar_generations(unsigned int _num_alveolar_generations) 
											{ num_alveolar_generations = _num_alveolar_generations; };
		unsigned int get_num_alveolar_generations() { return num_alveolar_generations; };

		void add_daughter(unsigned int daughter) { daughters.push_back(daughter); };
		std::vector<unsigned int> get_daughters() { return daughters; };
		int get_daughter_1();	
		bool has_daughter_1();	// same as has_daughters really

		void set_parent(int _parent) { parent = _parent; };
		int get_parent();	
		bool has_parent();


		void set_node_1(Point _node_1) { node_1 = _node_1; };
		Point get_node_1() { return node_1; };		

		void set_node_2(Point _node_2) { node_2 = _node_2; };
		Point get_node_2() { return node_2; };	

		void set_poiseuille(bool _poiseuille) { poiseuille = _poiseuille; };
		bool get_poiseuille() { return poiseuille; };	

		void set_local_elem_number(unsigned int _local_elem_number) { local_elem_number = _local_elem_number; };
		unsigned int get_local_elem_number() { return local_elem_number; };	

		void set_tree_number(unsigned int _tree_number) { tree_number = _tree_number; };
		unsigned int get_tree_number() { return tree_number; };	

		void set_flow(double _flow) { flow = _flow; };
		double get_flow() { return flow; };		

		void set_pressure_diff(double _pressure_diff) { pressure_diff = _pressure_diff; };
		double get_pressure_diff() { return pressure_diff; };		


		void set_velocity(double _velocity) { velocity = _velocity; };
		double get_velocity() { return velocity; };

		void set_previous_flow_rate(double _previous_flow_rate) { previous_flow_rate = _previous_flow_rate; };
		double get_previous_flow_rate() { return previous_flow_rate; };

		void set_previous_velocity(double _previous_velocity) { previous_velocity = _previous_velocity; };
		double get_previous_velocity() { return previous_velocity; };


		void set_velocity_2d(double _velocity_2d) { velocity_2d = _velocity_2d; };
		double get_velocity_2d() { return velocity_2d; };

		void set_previous_velocity_2d(double _previous_velocity_2d) { previous_velocity_2d = _previous_velocity_2d; };
		double get_previous_velocity_2d() { return previous_velocity_2d; };

		//move all element numbers by an amount
		void move_element_numbers(unsigned int amount);	

		void print();
		void print_concise();

		void add_particle_fraction(double amount, double entry_time);

		void remove_particle_fraction(unsigned int particle_fraction_number);

		unsigned int num_particle_fractions();

		double get_distance(unsigned int particle_fraction_number) { return particle_fraction_distance[particle_fraction_number]; };

		void set_distance(unsigned int particle_fraction_number, double distance) { particle_fraction_distance[particle_fraction_number] = distance; };

		double get_entry_time(unsigned int particle_fraction_number) { return particle_fraction_entry_time[particle_fraction_number]; };

		double get_particle_fraction(unsigned int particle_fraction_number) { return particle_fraction_amount[particle_fraction_number]; };

		double get_total_particle_fraction();
		void set_total_particle_fraction(double _particle_fraction);

		double get_cumulative_particle_fraction() { return cumulative_particle_fraction; };

		void set_deposition_probability_total(double _deposition_probability_total) { deposition_probability_total = _deposition_probability_total; };
		double get_deposition_probability_total() { return deposition_probability_total; };

		void set_deposition_probability_sed(double _deposition_probability_sed) { deposition_probability_sed = _deposition_probability_sed; };
		double get_deposition_probability_sed() { return deposition_probability_sed; };

		void set_deposition_probability_imp(double _deposition_probability_imp) { deposition_probability_imp = _deposition_probability_imp; };
		double get_deposition_probability_imp() { return deposition_probability_imp; };

		void set_deposition_probability_dif(double _deposition_probability_dif) { deposition_probability_dif = _deposition_probability_dif; };
		double get_deposition_probability_dif() { return deposition_probability_dif; };

		void set_stokes_number(double _stokes_number) { stokes_number = _stokes_number; };
		double get_stokes_number() { return stokes_number; };


		void set_deposition_fraction_total(double _deposition_fraction_total) { deposition_fraction_total = _deposition_fraction_total; };
		void add_deposition_fraction_total(double _deposition_fraction_total) { deposition_fraction_total += _deposition_fraction_total; };
		double get_deposition_fraction_total() { return deposition_fraction_total; };

		void set_deposition_fraction_sed(double _deposition_fraction_sed) { deposition_fraction_sed = _deposition_fraction_sed; };
		void add_deposition_fraction_sed(double _deposition_fraction_sed) { deposition_fraction_sed += _deposition_fraction_sed; };
		double get_deposition_fraction_sed() { return deposition_fraction_sed; };

		void set_deposition_fraction_imp(double _deposition_fraction_imp) { deposition_fraction_imp = _deposition_fraction_imp; };
		void add_deposition_fraction_imp(double _deposition_fraction_imp) { deposition_fraction_imp += _deposition_fraction_imp; };
		double get_deposition_fraction_imp() { return deposition_fraction_imp; };

		void set_deposition_fraction_dif(double _deposition_fraction_dif) { deposition_fraction_dif = _deposition_fraction_dif; };
		void add_deposition_fraction_dif(double _deposition_fraction_dif) { deposition_fraction_dif += _deposition_fraction_dif; };
		double get_deposition_fraction_dif() { return deposition_fraction_dif; };

		void set_terminal_exit_fraction(double _terminal_exit_fraction) { terminal_exit_fraction = _terminal_exit_fraction; };
		void add_terminal_exit_fraction(double _terminal_exit_fraction) { terminal_exit_fraction += _terminal_exit_fraction; };
		double get_terminal_exit_fraction() { return terminal_exit_fraction; };


		// symmetric airway data
		void set_is_symmetric_airway(bool _is_symmetric_airway) { is_symmetric_airway = _is_symmetric_airway; };
		bool get_is_symmetric_airway() { return is_symmetric_airway; };

		// num symmetric airways (so that can calculate how many particles etc)
		void set_num_symmetric_airways(unsigned int _num_symmetric_airways) { num_symmetric_airways = _num_symmetric_airways; };
		unsigned int get_num_symmetric_airways() { return num_symmetric_airways; };

		// bifurcation angle (e.g. useful for symmetric airways)
		// NOTE: only implmented for symmetric airways
		void set_branching_angle(double _branching_angle);
		double get_branching_angle();

		// gravity angle (e.g. useful for symmetric airways)
		// NOTE: only implmented for symmetric airways
		void set_gravity_angle(double _gravity_angle);
		double get_gravity_angle();

		// the acinar id (from 0)
		void set_acinar_id(int _acinar_id) { acinar_id = _acinar_id; };
		int get_acinar_id() { return acinar_id; };

		// the acinar elem number in the full mesh
		void set_acinar_elem(int _acinar_elem) { acinar_elem = _acinar_elem; };
		int get_acinar_elem() { return acinar_elem; };

		// the acinar elem number in the full mesh
		void set_acinar_output_elem(int _acinar_output_elem) { acinar_output_elem = _acinar_output_elem; };
		int get_acinar_output_elem() { return acinar_output_elem; };

		// the resistance of the airway
		void set_resistance(double _resistance) { resistance = _resistance; };
		double get_resistance() { return resistance; };

		// the poiseuille resistance of the airway
		void set_poiseuille_resistance(double _poiseuille_resistance) { poiseuille_resistance = _poiseuille_resistance; };
		double get_poiseuille_resistance() { return poiseuille_resistance; };

		// the wall thickness of the airway
		void set_wall_thickness(double _wall_thickness) { wall_thickness = _wall_thickness; };
		double get_wall_thickness() { return wall_thickness; };

		// the compliance of the airway (not the acinus)
		void set_airway_compliance(double _airway_compliance) { airway_compliance = _airway_compliance; };
		double get_airway_compliance() { return airway_compliance; };

		// the inertance
		void set_inertance(double _inertance) { inertance = _inertance; };
		double get_inertance() { return inertance; };

		// the local acinar compliance (if has an acinus, else -1)
		void set_local_acinar_compliance(double _local_acinar_compliance) { local_acinar_compliance = _local_acinar_compliance; };
		double get_local_acinar_compliance() { return local_acinar_compliance; };



	private:
		
		int generation;
		bool is_daughter_1;
		std::vector<unsigned int> siblings;
		int primary_sibling;
		double length;
		double radius;
		int order;
		double flow_rate;
		double efficiency;
		unsigned int num_alveolar_generations;
		std::vector<unsigned int> daughters;
		int parent;
		Point node_1;
		Point node_2;
		bool poiseuille;
		int local_elem_number;
		unsigned int tree_number;
		double flow;
		double pressure_diff;

		double previous_flow_rate;
		double velocity;
		double previous_velocity;
		double velocity_2d;
		double previous_velocity_2d;

		double deposition_probability_total;
		double deposition_probability_sed;
		double deposition_probability_imp;
		double deposition_probability_dif;

		double deposition_fraction_total;
		double deposition_fraction_sed;
		double deposition_fraction_imp;
		double deposition_fraction_dif;

		double particle_fraction;
		double cumulative_particle_fraction;

		double stokes_number;

		double terminal_exit_fraction;

		std::vector<double> particle_fraction_amount;
		std::vector<double> particle_fraction_entry_time;
		std::vector<double> particle_fraction_distance;

		bool is_symmetric_airway;
		unsigned int num_symmetric_airways;

		double branching_angle;
		double gravity_angle;

		int acinar_id;	// -1 if no acinus
		int acinar_elem;	// -1 if no acinus
		int acinar_output_elem;	// -1 if no acinus

		double resistance;	// updated after each assembly
		double poiseuille_resistance;	// set in the beginning
		double wall_thickness;	// set in the beginning
		double airway_compliance; // set in the beginning
		double inertance; // set in the beginning
		double local_acinar_compliance; // -1 if no acinus

};
#endif //__branch_h__

