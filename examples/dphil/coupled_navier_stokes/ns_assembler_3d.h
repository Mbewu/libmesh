
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

#include "surface_boundary.h"

#include "forcing_function.h"
#include "airway.h"

#include "exact_solution_velocity.h"
#include "exact_solution_pressure.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifndef __nsassembler3d_h__
#define __nsassembler3d_h__

// headers omitted for brevity
class NSAssembler3D : public System::Assembly
{
	public:
		NSAssembler3D (EquationSystems& es_in, std::vector<SurfaceBoundary* >& _surface_boundaries, std::vector<unsigned int> _subdomains_3d, unsigned int _n_initial_3d_elem) :
			es (&es_in),subdomains_3d(_subdomains_3d), n_initial_3d_elem(_n_initial_3d_elem),pressure_coupled(false),estimating_error(false),threed(true)
		
		{
			//rather insist on calling this yourself before each assembly so that you 
			//have the flexibility to change boundary conditions
			//init_bc ();  
			threed = es->parameters.get<bool>("threed");

			surface_boundaries = &_surface_boundaries;

		}

		// initialise bc information
		virtual void init_bc (std::vector<unsigned int> boundary_ids,
											std::vector<double> pressure_values,
											std::vector<double> previous_flux_values,
											std::vector<double> previous_previous_flux_values);

		virtual void assemble ()
		{
			ErrorVector error;
			this->assemble(error);
			if(es->parameters.get<unsigned int>("preconditioner_type"))
			{
				this->assemble_preconditioner();
			}
		}

		virtual void assemble (ErrorVector& error_vector) = 0;			//must be implemented by other class
		virtual void assemble_preconditioner ();								//must be implemented by other class (may not be though)
		virtual void calculate_peclet_number (ErrorVector& error_vector);
		virtual double calculate_flux (const int boundary_id);
		virtual double calculate_pressure (const int boundary_id);
		virtual double calculate_l2_norm(const unsigned int var_to_calc, bool reference=false);
		virtual void set_pressure_coupled (bool _pressure_coupled) { pressure_coupled = _pressure_coupled;};
		virtual void find_1d_boundary_nodes();	//here we set the primary_pressure_boundary_nodes_1d etc
		virtual void estimate_error(ErrorVector& _error_vector);	//here we set the primary_pressure_boundary_nodes_1d etc
		virtual void set_subdomains_1d (std::vector<unsigned int> _subdomains_1d) { subdomains_1d = _subdomains_1d;};
		virtual void set_airway_data (std::vector<Airway>& _airway_data) { airway_data = _airway_data;};
	
	protected:
		EquationSystems* es;
		// assuming they are initialised or whatevs
		std::map<const int,std::string> 			bc_type;
		std::map<const int,double> 			bc_value;
		std::map<const int,double> 			interp_flow_bc_value;
		std::map<const int,std::vector<Real> > boundary_centre;
		std::map<const int,Real> 							boundary_radius;
		bool pressure_coupled;
		std::vector<dof_id_type> primary_pressure_boundary_nodes_1d;
		std::vector<dof_id_type> secondary_pressure_boundary_nodes_1d;
		std::vector<dof_id_type> primary_flux_boundary_nodes_1d;
		std::vector<dof_id_type> secondary_flux_boundary_nodes_1d;
		bool estimating_error;
		bool threed;
		std::vector<SurfaceBoundary* >* surface_boundaries;
		std::vector<unsigned int> subdomains_3d;
		std::vector<unsigned int> subdomains_1d;
		int n_initial_3d_elem;

		std::vector<Airway> airway_data;

};

#endif //__nsassembler3d_assembly_h__

