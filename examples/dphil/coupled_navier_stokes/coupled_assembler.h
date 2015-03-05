
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

#include "ns_assembler_3d.h"
#include "navier_stokes_assembler.h"



// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifndef __coupled_assembler_h__
#define __coupled_assembler_h__

// headers omitted for brevity
class CoupledAssembler : public System::Assembly
{
	public:
		CoupledAssembler (NSAssembler3D& _ns_assembler_3d, NavierStokesAssembler& _ns_assembler) :
			p_ns_assembler_3d(&_ns_assembler_3d), p_ns_assembler(&_ns_assembler)
		{
 
		}

		void assemble() 
		{
			p_ns_assembler_3d->assemble();
			p_ns_assembler->assemble();
		}
	
	private:
		NSAssembler3D* p_ns_assembler_3d;
		NavierStokesAssembler* p_ns_assembler;
		
};

#endif //__picard_assembly_h__

