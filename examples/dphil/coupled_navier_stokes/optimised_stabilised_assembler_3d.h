
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

#include "libmesh/error_vector.h"

#include "ns_assembler_3d.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifndef __optimisedstabilisedassembler3d_h__
#define __optimisedstabilisedassembler3d_h__

// headers omitted for brevity
class OptimisedStabilisedAssembler3D : public NSAssembler3D
{
	public:
		OptimisedStabilisedAssembler3D (EquationSystems& es_in, 
																	std::vector<SurfaceBoundary* >& _surface_boundaries, 
																	std::vector<unsigned int> subdomains_3d, 
																	//unsigned int n_initial_3d_elem) :
																	std::vector<unsigned int> elem_to_airway) :
			NSAssembler3D (es_in,_surface_boundaries,subdomains_3d, elem_to_airway)
		{}

		using NSAssembler3D::assemble;
		virtual void assemble (ErrorVector& error_vector);

};

#endif //__optimisedstabilisedassembler3d_h__

