
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

#ifndef __particles3d_h__
#define __particles3d_h__

// headers omitted for brevity
class Particle
{
	public:
		Particle (EquationSystems& es_in,Point& p,const Elem* element, unsigned int _particle_id, double _perturbation_magnitude = 1e-5);

		// initialise particles
		Point get_position (){ return position;};

		// move particles over one timestep
		const Elem* get_current_elem ()
		{ 
			if(!exited)
				return current_elem;
			else
			{
				std::cout << "ERROR: particle has left domain and therefore has no current elem... EXITING" << std::endl;
				std::exit(0);
			}
		};

		// move particles over one timestep
		int get_exit_surface () { return exit_surface; };

		// move particles over one timestep
		bool is_on_wall() { return on_wall; };

		// move particles over one timestep
		bool has_exited() { return exited; };

		// move particles over one timestep
		double get_exit_time() { return time_exited;};

		// move particles over one timestep
		double get_entrance_time() { return entrance_time;};

		// move particles over one timestep
		double get_particle_id() { return particle_id;};

		double get_maximum_timestep (){ return maximum_timestep; };

		NumberVectorValue get_current_velocity() { return current_velocity; };

		void try_and_move();

		int get_status()
		{
			if(on_wall)
				return 1;
			else if(exited)
				return 2;
			else if(broken)
				return 3;
			else
				return 0;
		};
	
	private:

		NumberVectorValue compute_velocity ();

		NumberVectorValue compute_particle_velocity (NumberVectorValue velocity);

		bool goes_through_side(Point new_position,int side_number,const Elem* element, double old_s_param, double& s_param);

		void check_neighbors(const Elem* element,Point& new_position, double old_s_param, std::vector<unsigned int>& elements_checked, bool& element_found);

		void move();

		double particle_reynolds_number ();

		double drag_coeff (double reynolds_number);

		
		
		EquationSystems* es;
		Point position;
		const Elem* current_elem;		//-1 if outside of domain.
		bool on_wall;		//true if on a wall
		bool exited;		//true if exited domain
		int exit_surface;		//surface through which particle passed
		double maximum_timestep;		//this is the maximum timestep allowed. (equiv to moving an element length)
		NumberVectorValue current_velocity;
		NumberVectorValue current_particle_velocity;
		double time_exited;
		bool threed;
		unsigned int particle_id;
		double entrance_time;
		bool broken;


		double perturbation_magnitude;
		double constant_drag_force;
		double cunningham_correction_factor;

		bool sedimentation;
		bool impaction;
		bool drag;


};

#endif //__particles3d_h__

