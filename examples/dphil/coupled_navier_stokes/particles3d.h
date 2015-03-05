
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
#include "particle.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifndef __particles3d_h__
#define __particles3d_h__

// headers omitted for brevity
class Particles3D
{
	public:
		Particles3D (EquationSystems& es_in, std::vector<SurfaceBoundary* >& _surface_boundaries) :
			es (&es_in),threed(true)
		
		{
			//rather insist on calling this yourself before each assembly so that you 
			//have the flexibility to change boundary conditions
			//init_bc ();  
			threed = es->parameters.get<bool>("threed");

			surface_boundaries = &_surface_boundaries;

		}

		// initialise particles
		void init_particles ();

		// move particles over one timestep
		void move_particles ();

	
	protected:
		EquationSystems* es;
		bool threed;
		std::vector<SurfaceBoundary* >* surface_boundaries;
		std::vector<Particle> particles;

};

#endif //__particles3d_h__

