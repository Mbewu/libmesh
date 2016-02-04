
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

#include "airway.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifndef __navier_stokes_assembler_h__
#define __navier_stokes_assembler_h__

// headers omitted for brevity
class NavierStokesAssembler : public System::Assembly
{
	public:
		NavierStokesAssembler (EquationSystems& es_in, std::vector<Airway>& _airway_data,	std::vector<unsigned int> _subdomains_1d, unsigned int _n_initial_3d_elem) :
			es (&es_in), subdomains_1d(_subdomains_1d), airway_data(_airway_data), n_initial_3d_elem(_n_initial_3d_elem), coupled(false)
		{
		}

		// initialise bc information
		void assemble ();
		void assemble_stokes_steady_0D ();
		void assemble_stokes_1D ();
		double calculate_flux (const int boundary_id, const int mid_mesh_element=-1);
		double calculate_pressure (const int boundary_id, const int mid_mesh_element=-1);
		void init_bc(std::vector<double> _flux_values, std::vector<double> _pressure_values = std::vector<double>());
		void set_coupled(bool _coupled) { coupled = _coupled; };
		void add_to_matrices(TransientLinearImplicitSystem * system, unsigned int row_number, unsigned int col_number, double value);
	
	private:
		EquationSystems* es;
		std::vector<Airway>& airway_data;
		int n_initial_3d_elem;
		std::vector<double> pressure_values;
		std::vector<double> flux_values;
		bool coupled;
		std::vector<unsigned int> subdomains_1d;
		bool preconditioner;
};

#endif //__navier_stokes_assembler_assembly_h__

