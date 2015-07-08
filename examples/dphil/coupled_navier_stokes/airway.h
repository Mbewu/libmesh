
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

		void set_order(int _order) { order = _order; };
		int get_order();

		void set_flow_rate(double _flow_rate) { flow_rate = _flow_rate; };
		double get_flow_rate() { return flow_rate; };

		void set_efficiency(double _efficiency) { efficiency = _efficiency; };
		double get_efficiency() { return efficiency; };

		void set_num_alveolar_generations(unsigned int _num_alveolar_generations) 
											{ num_alveolar_generations = _num_alveolar_generations; };
		unsigned int get_num_alveolar_generations() { return num_alveolar_generations; };

		void add_daughter(unsigned int daughter) { daughters.push_back(daughter); };
		std::vector<unsigned int> get_daughters() { return daughters; };
		int get_daughter_1();	
		bool has_daughter_1();	

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

		//move all element numbers by an amount
		void move_element_numbers(unsigned int amount);	

		void print();
		void print_concise();


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
		unsigned int local_elem_number;
		unsigned int tree_number;
		double flow;
		double pressure_diff;

};
#endif //__branch_h__

