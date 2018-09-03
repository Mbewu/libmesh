/* The libMesh Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

 // <h1>Systems Example 2 - Unsteady Nonlinear Navier-Stokes</h1>
 //
 // This example shows how a simple, unsteady, nonlinear system of equations
 // can be solved in parallel.  The system of equations are the familiar
 // Navier-Stokes equations for low-speed incompressible fluid flow.  This
 // example introduces the concept of the inner nonlinear loop for each
 // timestep, and requires a good deal of linear algebra number-crunching
 // at each step.  If you have a ExodusII viewer such as ParaView installed,
 // the script movie.sh in this directory will also take appropriate screen
 // shots of each of the solution files in the time sequence.  These rgb files
 // can then be animated with the "animate" utility of ImageMagick if it is
 // installed on your system.  On a PIII 1GHz machine in debug mode, this
 // example takes a little over a minute to run.  If you would like to see
 // a more detailed time history, or compute more timesteps, that is certainly
 // possible by changing the n_timesteps and dt variables below.


// C++ include files that we need
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <math.h>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "exodusII_io_extended.h"
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
#include "libmesh/zero_function.h"
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"


// local includes
#include "augment_sparsity_on_interface.h"

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

// Defines the MeshData class, which allows you to store
// data about the mesh when reading in files, etc.
#include "libmesh/mesh_data.h"

#include "libmesh/gmv_io.h"
#include "libmesh/serial_mesh.h"

#include "picard.h"
#include "optimised_stabilised_assembler_3d.h"
#include "ns_assembler_3d.h"

// for fieldsplit stuff
#include "libmesh/petscdmlibmesh.h"
//#include "petsc_dm_nonlinear_solver.C"

#include "libmesh/petsc_linear_solver.h"

//mesh refinement
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"
#include "libmesh/uniform_refinement_estimator.h"
#include "z_position_flagging.h"
#include "xyz_position_flagging.h"
#include "ball_boundary_flagging.h"


//for mesh generation
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"

//1d assembler
#include "navier_stokes_assembler.h"

//coupled assembler
#include "coupled_assembler.h"

#include "augment_sparsity_on_interface.h"
//#include "enriched_analytic_function.h"
#include "inflow_boundary_function.h"
#include "moving_pipe_boundary_function.h"
#include "lid_driven_cavity_boundary_function.h"

//for making variables active only on a subdomain
#include "libmesh/system_subset_by_subdomain.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/mesh_modification.h"

//for file system stuff
#include <dirent.h>

#include "libmesh/fem_function_base.h"

#include "surface_boundary.h"

#include "exact_solution_velocity.h"
#include "exact_solution_pressure.h"
#include "airway.h"
#include "libmesh/exact_solution.h"

#include "particle.h"
#include "hofmann_particle_deposition.h"
#include "james_particle_deposition.h"


// Needed for cast to petsc matrix
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"



#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_base.h"

// for element type names
#include "libmesh/elem_type.h"
#include "libmesh/string_to_enum.h"

#include "augment_sparsity_moghadam.h"

#include "augment_sparsity_acinar.h"
#include "libmesh/node_elem.h"

#include "augment_sparsity_combined.h"

//#include <boost/bind.hpp>

// *** stuff for shell pc


typedef struct {
	Vec inner_solution;
} InitialGuessCtx;

/* Define context for user-provided preconditioner */
typedef struct {
  Mat pressure_mass_matrix;
  Mat pressure_laplacian_matrix;
  Mat pressure_laplacian_preconditioner;
  Mat pressure_convection_diffusion_matrix;
  Mat velocity_matrix;
  Mat S_approx;
  KSP inner_mass_ksp;
  KSP inner_lap_ksp;
  KSP inner_schur_ksp;
  KSP inner_velocity_ksp;
  KSP outer_ksp;
  Vec temp_vec;
  Vec temp_vec_2;
  Vec temp_vec_3;
  Vec lsc_scale;
  Mat lsc_laplacian_matrix;
  Vec lsc_stab_alpha_D_inv;
	bool neg_schur;
  PetscInt total_iterations;
  Mat Bt;
  InitialGuessCtx *initial_guess_ctx;
} NSShellPC;

typedef struct {
  Mat pressure_mass_matrix;
  Mat pressure_laplacian_matrix;
  Mat pressure_laplacian_preconditioner;
  Mat pressure_convection_diffusion_matrix;
  Mat velocity_matrix;
  Mat S_approx;
  KSP inner_mass_ksp;
  KSP inner_lap_ksp;
  KSP inner_schur_ksp;
  KSP inner_velocity_ksp;
  KSP outer_ksp;
  Vec temp_vec;
  Vec temp_vec_2;
  Vec temp_vec_3;
  Vec lsc_scale;
  Mat lsc_laplacian_matrix;
  Vec lsc_stab_alpha_D_inv;
	bool neg_schur;
  PetscInt total_iterations;
  Mat Bt;
  InitialGuessCtx *initial_guess_ctx;
} MonoShellPC;

typedef struct {
  Mat velocity_matrix;
  Mat S_approx;
  Mat A;
  Mat Bt;
  Mat B;
  Mat C;
  IS velocity_is;
  IS pressure_is;
  KSP inner_schur_ksp;
  KSP inner_velocity_ksp;
  Vec temp_vec;
  Vec temp_vec_2;
  Vec temp_vec_3;
} SIMPLEShellPC;

/* Define context for user-provided matrix A_p = BQ^-1Bt */

typedef struct {
  //Mat velocity_diag;
  Mat velocity_mass_matrix;
  Mat pressure_laplacian_matrix;
  //Mat b_matrix;
  //Mat bt_matrix;
	KSP schur_ksp;
	PetscInt m;
	PetscInt n;
} PCD2ShellMatrixCtx;



/* Define context for user-provided matrix S_tilde = C - BHBt */

typedef struct {
  //Mat velocity_diag;
  Mat moghadam_velocity_preconditioner_matrix;
  Mat moghadam_b_matrix;
  Mat moghadam_bt_matrix;
  Mat moghadam_c_matrix;
  //Mat b_matrix;
  //Mat bt_matrix;
	//KSP schur_ksp;
	//PetscInt m;
	//PetscInt n;
} MoghadamSchurShellMatrixCtx;

/* Define context for user-provided preconditioner */
typedef struct {
  Mat moghadam_velocity_preconditioner_matrix;
} MoghadamVelocityPCCtx;

typedef struct {
	PetscInt total_velocity_iterations;
	PetscInt total_convection_diffusion_iterations;
} MonolithicMonitorCtx;





/* Declare routines for user-provided preconditioner */
extern PetscErrorCode ShellPCCreate(NSShellPC**);
extern PetscErrorCode ShellPCCreate(MonoShellPC**);
extern PetscErrorCode ShellPCCreate(SIMPLEShellPC**);
extern PetscErrorCode PressureShellPCSetUp(PC,Mat,KSP);
extern PetscErrorCode PressureShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode PCDShellPCSetUp(PC,Mat,Mat,Mat,Mat,KSP);
extern PetscErrorCode PCDShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode PCD2ShellPCSetUp(PC,Mat,Mat,Mat,Mat,KSP);
extern PetscErrorCode PCD2ShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode MonolithicShellPCSetUp(PC,Mat,KSP,KSP);
extern PetscErrorCode Monolithic3ShellPCSetUp(PC,Mat,KSP,KSP,bool,bool,double);
extern PetscErrorCode Monolithic4ShellPCSetUp(PC,Mat,Vec,Vec,KSP,KSP,bool,bool,double);
extern PetscErrorCode Monolithic2ShellPCSetUp(PC,Mat,Vec,Vec,KSP,KSP);
extern PetscErrorCode MonolithicShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode Monolithic2ShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode LSCShellPCSetUp(PC,KSP);
extern PetscErrorCode LSCShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode LSCScaledShellPCSetUp(PC,Mat,KSP);
extern PetscErrorCode LSCScaledShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode LSCScaledStabilisedShellPCSetUp(PC,Mat,KSP,bool);
extern PetscErrorCode LSCScaledStabilisedShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode SIMPLEShellPCSetUp(PC,IS,IS,KSP);
extern PetscErrorCode SIMPLECShellPCSetUp(PC,IS,IS,KSP);
extern PetscErrorCode SIMPLEShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode SIMPLERShellPCApply(PC,Vec x,Vec y);
extern PetscErrorCode NSShellDestroy(PC);
extern PetscErrorCode MonoShellDestroy(PC);
extern PetscErrorCode SIMPLEShellDestroy(PC);

extern PetscErrorCode custom_outer_monitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy);
extern PetscErrorCode custom_inner_pressure_monitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy);
extern PetscErrorCode custom_inner_velocity_monitor(KSP ksp, PetscInt n, PetscReal rnorm, void *dummy);
extern PetscErrorCode compute_initial_guess(KSP ksp, Vec x, void *ctx);


extern PetscErrorCode compute_matrix(KSP,Mat,Mat,void*);
extern PetscErrorCode compute_rhs(KSP,Vec,void*);


/* Declare routines for user-provided A_p = BQ^-1Bt */
extern PetscErrorCode MatShellMultFull(Mat,Vec,Vec);
//extern PetscErrorCode MatShellMultDiag(Mat,Vec,Vec);
extern PetscErrorCode MoghadamSchurMatShellMultFull(Mat,Vec,Vec);
extern PetscErrorCode MoghadamNoPCSchurMatShellMultFull(Mat,Vec,Vec);
extern PetscErrorCode MoghadamExactPCSchurMatShellMultFull(Mat,Vec,Vec);
extern PetscErrorCode MoghadamVelocityPCApply(PC,Vec,Vec);
extern PetscErrorCode ShellPCCreate(MoghadamVelocityPCCtx**);
extern PetscErrorCode MoghadamVelocityPCDestroy(PC);
extern PetscErrorCode MoghadamVelocityPCSetUp(PC,Mat);

// Bring in everything from the libMesh namespace
using namespace libMesh;





// ********************** CLASS DEFINITION ****************************** //

class NavierStokesCoupled
{

public:
	//constructor for the class, just used to store 
	//local variables between functions etc
	NavierStokesCoupled(LibMeshInit & init, std::string input_file, std::string input_file_particle, GetPot& command_line);

	// read in parameters from file and give them to the parameters object.
	int read_parameters();

	void output_parameters();

	void print_parameters();

	void output_command_line_options();

	void set_auto_fieldsplit_parameters();

	void read_particle_parameters();

	void output_particle_parameters();

	void setup_read_timesteps();

	// output_mesh true means we don't set up extra crap
	void setup_3d_mesh(EquationSystems* _es, Mesh& _mesh, bool output_mesh=false);

	void setup_3d_system(System * system, bool output_system=false);

	void prerefine_3d_mesh();

	void setup_1d_mesh ();

	void read_1d_mesh (bool centerlines_mesh=false);

	void generate_1d_mesh ();

	void setup_1d_system(System * system, bool output_system=false, bool centrelines=false);

	void setup_acinar_system(System * system, bool output_system=false, bool centrelines=false);

	void output_sim_data(bool header=false);

	void output_linear_iteration_count(bool header=false);

	void output_3d_particle_data();

	void output_particle_data_old(bool header=false);

	void calculate_3d_particle_data();

	void calculate_0d_particle_data();

	void output_deposition_data(bool header=false);

	bool no_particles_in_motion();

	void print_flux_and_pressure();

	// write the solution at a specific time and time_step, 
	// backup=TRUE if backup is to be written
	void write_3d_solution();

	void write_3d_particles();

	void write_1d_solution();

	void write_acinar_solution();

	void write_efficiency_solution();

	void update_times();

	void update_time_scaling();

	int solve_1d_system(TransientLinearImplicitSystem * system);

	void set_radii(System * system, bool centrelines=false);

	void set_poiseuille();

	void set_efficiency();

	void calculate_1d_boundary_values(bool mono_0d_only=false);

	//returns true if wanna break iteration
	bool solve_3d_system_iteration(TransientLinearImplicitSystem * system);

	void assemble_3d_system(TransientLinearImplicitSystem * system);

	double solve_3d_system(TransientLinearImplicitSystem * system);

	void solve_3d_clean_up(TransientLinearImplicitSystem * system, Mat& B);

	double solve_3d_system_moghadam(TransientLinearImplicitSystem * system);

	void solve_3d_clean_up_moghadam_matrices(TransientLinearImplicitSystem * system);

	void solve_3d_clean_up_moghadam_vectors(TransientLinearImplicitSystem * system);

	double solve_3d_system_residual(TransientLinearImplicitSystem * system);

	void calculate_3d_boundary_values();

	void adaptively_refine();

	void end_timestep_admin();

	unsigned int add_1d_tree_to_mesh(Mesh& _mesh, std::vector<Point>& vertices, std::vector<std::vector<unsigned int> >& cell_vertices,unsigned int subdomain_id, unsigned int boundary_id, bool full_mesh = false);

	void add_acinar_to_mesh(Mesh& _mesh, std::vector<Point>& vertices, std::vector<std::vector<unsigned int> >& cell_vertices,unsigned int subdomain_id, unsigned int airway_elem_id_start, bool full_mesh = false);

	void create_1d_tree(std::vector<Point>& vertices, std::vector<std::vector<unsigned int> >& cell_vertices, unsigned int num_generations);

	void calculate_num_alveloar_generations_for_tree();

	// convert the 1d part of a vector from monomial to nodal
	void convert_1d_solution_monomial_to_nodal();

	// convert the 1d part of a vector from nodal to monomial
	void convert_1d_solution_nodal_to_monomial();

	void plot_error(ExodusII_IO_Extended& io);

	void write_elem_pid_1d(ExodusII_IO_Extended& io);

	void write_elem_pid_3d(ExodusII_IO_Extended& io);

	void set_elem_proc_id_3d();
	void set_elem_proc_id_1d(System * system, bool centrelines=false);

	double set_double_parameter(GetPot _infile, std::string name, double default_value);
	unsigned int set_unsigned_int_parameter(GetPot _infile, std::string name, unsigned int default_value);
	int set_int_parameter(GetPot _infile, std::string name, int default_value);
	bool set_bool_parameter(GetPot _infile, std::string name, bool default_value);
	std::string set_string_parameter(GetPot _infile, std::string name, std::string default_value);
	std::string set_string_command_line_parameter(std::string name, std::string default_value);

	double scaled_time();

	void scale_3d_solution_vector(double velocity_scaling,double pressure_scaling);

	void scale_1d_solution_vector(double flux_scaling,double pressure_scaling);

	void init_dof_variable_vector(System * system, std::vector<int>& _dof_variable_type);

	void init_dof_variable_vectors();

	//void setup_results_vector(std::string results_folder, 
	//			std::vector<NumericVector<Number>* >& results_vector);		//setup vectors that stores the results from two previous simulations

	void output_poiseuille_resistance_per_generation();
	//okay want some functions that output scaled quantities
	//we can do this by constructing temporary explicit equation systems?

	void init_particles();

	double init_0d_particles();

	void move_particles();


	void setup_variable_scalings_3D();

	void setup_variable_scalings_1D();

	void setup_variable_scalings_acinar();

	void output_coupling_points();



	void construct_petsc_options_string();

	int setup_preconditioners(TransientLinearImplicitSystem * system, Mat& B);

	int libmesh_jacobi_pre(TransientLinearImplicitSystem * system);

	int libmesh_jacobi_post(TransientLinearImplicitSystem * system);

	int setup_moghadam_matrices(TransientLinearImplicitSystem * system);

	int setup_moghadam_solvers(TransientLinearImplicitSystem * system);

	int moghadam_libmesh_jacobi_pre(TransientLinearImplicitSystem * system);

	int moghadam_libmesh_jacobi_pre_rhs(TransientLinearImplicitSystem * system);

	int moghadam_libmesh_jacobi_post(TransientLinearImplicitSystem * system);

	int setup_moghadam_velocity_1(TransientLinearImplicitSystem * system);

	int setup_moghadam_schur(TransientLinearImplicitSystem * system);

	int setup_moghadam_velocity_2(TransientLinearImplicitSystem * system);

	int update_moghadam_solution_vector(TransientLinearImplicitSystem * system);

	int check_moghadam_preconditioned_solution(TransientLinearImplicitSystem * system);

	int compute_and_output_eigenvalues(TransientLinearImplicitSystem * system);

	int setup_is_simple (TransientLinearImplicitSystem * sys, PC my_pc);

	int construct_schur_stokes_matrix(TransientLinearImplicitSystem * sys, KSP schur_ksp);

	void set_reynolds_number();

	void set_characteristic_length_0d();

	void set_characteristic_length_3d();


	int test_post_solve(TransientLinearImplicitSystem * sys);

	// read the input boundary conditions for an uncoupled simulation
	int read_input_boundary_conditions_3d0d();

	// read the input boundary conditions for a pressure simulation
	int read_input_boundary_conditions_3d();

	void calculate_1d_pressure_deriv_values();

	void update_3d_dirichlet_boundary_conditions();

	void petsc_clean_up();

	void output_logging();

	void read_old_solution_from_file(EquationSystems* _es, std::ostringstream& file_name, unsigned int read_time_step);
	
	void copy_3d_solution_output_to_program();

	void copy_3d_solution_program_to_output();

	void copy_1d_solution_output_to_program();

	void copy_1d_solution_program_to_output();

	void copy_acinar_solution_program_to_output();

	void read_old_solutions(unsigned int read_time_step, double read_time);

	void copy_back_1d_solutions();

	void calculate_local_particle_1d_flow_rate();

	void copy_back_3d_solutions();

	void calculate_local_particle_3d_flow_rate();

	// set the deposition fraction in the 2D/3D airways.
	void set_centrelines_deposition_data ();

	void write_centrelines ();

	void output_3d_airway_deposition_data ();

	void output_0d_airway_deposition_data ();

	void output_global_deposition_metrics (bool header=false);

	// this outputs the tree structure of the airways
	void output_airway_tree ();

	std::vector<int> calculate_3d_to_0d_airway_ids();

	void setup_pressure_boundary_conditions();

	double calculate_total_acinar_volume();

	void init_acinar_volumes();

	double calculate_total_terminal_area();

	double calculate_average_terminal_pressure();

	void set_resistance();

	void init_1d_airway_parameters();

	std::string get_converged_reason_string(int);

	void output_timings_to_file();

	void calculate_total_pressure();

private:

	unsigned int sim_type;	//0 - 3d, 1 - 1d, 2 - uncoupled, 3 - expl coupled
	std::vector<std::vector<double> > element_data;
	std::vector<Airway> airway_data;
	std::vector<std::vector<Airway> > terminal_airway_data;
	std::vector<std::vector<double> > centreline_element_data;	//element data of centrelines
	std::vector<Airway> centreline_airway_data;
	Mesh mesh;
	Mesh mesh_3d;	// 3d only mesh for coupled output
	Mesh mesh_1d;	// 0d only mesh for coupled output
	Mesh mesh_1d_terminal;	// 0d only mesh for coupled output
	Mesh mesh_acinar;	// 0d only mesh for coupled output
	Mesh mesh_1d_acinar;	// 0d only mesh for coupled output
	Mesh mesh_centrelines;	// 0d only mesh for centrelines output
	MeshData mesh_data;
	MeshRefinement mesh_refinement;
	MeshRefinement mesh_refinement_3d;
	MeshRefinement mesh_refinement_1d;
	AutoPtr<EquationSystems> es;
	AutoPtr<EquationSystems> es_3d;	// eq for 3d output
	AutoPtr<EquationSystems> es_1d; // eq for 0d output
	AutoPtr<EquationSystems> es_1d_terminal; // eq for 0d terminal output
	AutoPtr<EquationSystems> es_centrelines; // eq for centrelines output
	AutoPtr<EquationSystems> es_acinar; // eq for 0d output
	TransientLinearImplicitSystem * system_3d;
	TransientLinearImplicitSystem * system_1d;
	TransientExplicitSystem * system_3d_output;	// system for 3d output
	TransientExplicitSystem * system_1d_output;	// system for 0d output
	TransientExplicitSystem * system_acinar_output;	// system for 0d output
	TransientLinearImplicitSystem * system_coupled;
	TransientLinearImplicitSystem * system_neumann;
	TransientExplicitSystem * extra_1d_data_system;
	TransientExplicitSystem * extra_3d_data_system;
	TransientExplicitSystem * particle_deposition_system_1d;
	TransientExplicitSystem * system_centrelines;
	unsigned int t_step;	
	double time;
	double dt;
	double previous_dt;
	AutoPtr<NSAssembler3D> picard;	//instead of being picard, this should be 3d assembly function
	AutoPtr<NavierStokesAssembler> ns_assembler;
	AutoPtr<NavierStokesAssembler> ns_assembler_mono;
	AutoPtr<CoupledAssembler> coupled_assembler;
	AutoPtr<ErrorEstimator> error_estimator;
	ErrorVector error_vector;
	std::vector<double> re_vec;
	unsigned int reduce_dt;	// 0 - don't, 1 - because nonlin, 2 - because linear, 3 - because residual
	bool increase_dt;
	bool refine_mesh;
	unsigned int steps_since_last_dt_change;
	unsigned int unsteady;
	bool restart;
	PerfLog perf_log;
	PerfLog perf_log_move;
	AutoPtr<NumericVector<Number> > last_nonlinear_soln;	//could be pointer but too much admin;
	AutoPtr<NumericVector<Number> > last_1d_nonlinear_soln;	//could be pointer but too much admin
	AutoPtr<NumericVector<Number> > old_global_solution;	//could be pointer but too much admin
	AutoPtr<NumericVector<Number> > old_global_solution_1d;	//could be pointer but too much admin
	Real previous_nonlinear_residual;
	Real current_nonlinear_residual;
	Real current_1d_nonlinear_residual;
	bool sim_3d;
	bool sim_1d;
	std::ofstream output_file;
	std::ofstream eigenvalues_file;
	std::ofstream eigenvalues_approx_file;
	std::ofstream parameter_output_file;
	std::ofstream particle_output_file;
	std::ofstream deposition_output_file;
	std::ofstream linear_iterations_output_file;
	std::ofstream global_deposition_metrics_file;
	std::vector<unsigned int> boundary_ids;
	// these are not boundary conditions but values on the respective domains
	// calculated after the fact.
	// e.g. may want to apply 3d values as boundary conditions
	std::vector<double> pressure_values_3d;
	std::vector<double> flux_values_3d;
	std::vector<double> pressure_values_1d;
	std::vector<double> flux_values_1d;
	std::vector<double> previous_pressure_values_3d;
	std::vector<double> previous_previous_pressure_values_3d;
	std::vector<double> previous_flux_values_3d;
	std::vector<double> previous_previous_flux_values_3d;
	std::vector<double> previous_pressure_values_1d;
	std::vector<double> previous_flux_values_1d;
	std::vector<double> input_pressure_values_3d;
	std::vector<double> input_flux_values_1d;
	std::vector<double> input_flux_values_3d;
	std::vector<double> pressure_deriv_values_1d;
	std::vector<double> pressure_zero_flow_values_1d;
	std::vector<double> previous_dynamic_pressure_values_3d;
	std::vector<double> total_pressure_values_3d;
	AutoPtr<AugmentSparsityOnInterface> augment_sparsity;
	AutoPtr<AugmentSparsityMoghadam> augment_sparsity_moghadam;
	AutoPtr<AugmentSparsityAcinar> augment_sparsity_acinar;
	AutoPtr<AugmentSparsityCombined> augment_sparsity_combined_1d;
	AutoPtr<AugmentSparsityCombined> augment_sparsity_combined_3d;
	AutoPtr<AugmentSparsityCombined> augment_sparsity_combined_coupled;
	//std::string _input_file;
	std::string input_file;
	std::string input_file_particle;
	std::ostringstream output_folder;
	std::ostringstream restart_folder;
	GetPot infile;
	GetPot infileparticle;
	GetPot comm_line;
	bool exit_program;							//basically quits all the loops and such
	double time_scale_factor;
	std::vector<int> dof_variable_type_3d;	//for each dof, has appropriate variable id
	std::vector<int> dof_variable_type_1d;	//for each dof, has appropriate variable id
	std::vector<int> dof_variable_type_coupled;	//for each dof, has appropriate variable id
	std::vector<int> dof_variable_type_acinar_output;	//for each dof, has appropriate variable id
	std::vector<int> dof_variable_type_3d_output;	//for each dof, has appropriate variable id
	std::vector<int> dof_variable_type_1d_output;	//for each dof, has appropriate variable id
	std::vector<int> dof_variable_type_3d_extra;	//for each dof, has appropriate variable id
	//std::vector<NumericVector<Number>* > results_vector_1;		//could be a vector of pointers but too much admin
	//std::vector<NumericVector<Number>* > results_vector_2;		//could be a vector of pointers but too much admin
	bool threed;
	std::vector<unsigned int> subdomains_1d;		//list of subdomain ids that are 1d (i.e. 1d airways)
	std::vector<unsigned int> subdomains_1d_terminal;		//list of subdomain ids that are 1d (i.e. 1d airways)
	std::vector<unsigned int> subdomains_acinar;		//list of subdomain ids that are 1d (i.e. 1d airways)
	std::vector<unsigned int> subdomains_2d;		//list of subdomain ids that are 2d (i.e. neumann boundary)
	std::vector<unsigned int> subdomains_3d;		//list of subdomain ids that are 3d (i.e. 3d airways) 
	unsigned int total_nonlinear_iterations;		//list of subdomain ids that are 3d (i.e. 3d airways)
	unsigned int local_linear_iterations;		//list of subdomain ids that are 3d (i.e. 3d airways) 
	int total_gmres_linear_iterations;
	int total_cg_linear_iterations;
	int total_linear_iterations;
	int total_max_iterations;
	int total_max_gmres_iterations;
	int total_max_cg_iterations;
	int local_max_iterations;
	int local_max_residual_iterations;
	//SurfaceBoundary inflow_surface_boundary_object;
	std::vector<SurfaceBoundary* > surface_boundaries;
	unsigned int nonlinear_iteration;
	unsigned int nonlinear_iteration_1d;
	unsigned int num_linear_iterations;
	//edge num is idx, 0-elem_num 1-node1 2-node2 3-generation 4-order
	std::vector<std::vector<double> > edge_1d_data;
	unsigned int num_edges_1d;
	//node num is idx, 0-xcoord 1-ycoord 2-zcoord 3-radius
	std::vector<std::vector<double> > node_1d_data;
	unsigned int num_nodes_1d;
	std::vector<unsigned int> num_generations; //num_generations for each tree
	bool ic_set;
	unsigned int residual_linear_iteration;

	// acinar stuff
	std::vector<unsigned int> acinar_to_airway;	// acinar no to airway elem no
	std::vector<unsigned int> elem_to_acinar;	// element number (in full mesh) to acinar no
	std::vector<unsigned int> elem_to_airway;	// element number (in full mesh) to airway no
	unsigned int num_airways_in_full_mesh;

	// particle deposition stuff
	unsigned int particle_deposition;
	std::vector<double> timestep_sizes;	//vector of timestep sizes, 1st timestep is 0th element obv
	std::vector<Particle> particles_3D;

	std::vector<double> var_scalings_3D;
	std::vector<double> var_scalings_1D;
	std::vector<double> var_scalings_acinar;
	std::vector<double> total_efficiency;


	std::vector<int> subtree_starting_elements;	//ending elements of the 3d mesh and the subtree that they go to
	std::vector<int> boundary_id_to_tree_id;
	std::vector<int> tree_id_to_boundary_id;
	std::vector<double> coupling_pressure_values_1d;
	std::vector<double> coupling_flux_values_1d;

	std::vector<int> centreline_terminal_id_to_tree_id;
	std::vector<Point> centreline_points_1;
	std::vector<Point> centreline_points_2;

	unsigned int max_generations_3d;
	unsigned int max_generations_0d;
	unsigned int max_generations;

	// this is a vector of the coupling points for each tree
	std::vector<Point> coupling_points;

	std::vector<std::vector<double> > alveolar_efficiency_per_generation;	// [timestep][generation] alveolar efficiency of deposition summed over all branches in this generation
	std::vector<std::vector<double> > tb_efficiency_per_generation;	// [timestep][generation] tracheo-bronchial efficiency of deposition summed over all branches in this generation
	int total_particles_inhaled;
	int stokes_gmres_iterations;
	bool shell_pc_created;
	bool schur_shell_pc_created;
	bool mono_shell_pc_created;
	bool first_3d_write;
	bool first_1d_write;
	bool init_names_done;

	std::vector<std::vector<double> > all_input_pressure_values_3d;
	std::vector<std::vector<double> > all_input_flux_values_0d;
	std::vector<std::vector<double> > all_input_flux_values_3d;
	std::string input_boundary_conditions_filename;

	// james particle deposition stuff
	JamesParticleDeposition james_particle_deposition_object;
	
	std::vector<std::vector<unsigned int>> airway_elem_id_starts;
	std::vector<unsigned int> terminal_airway_elem_id_starts;	// only need one number cause symmetric
	std::vector<double> total_fraction_exited_3d_surface;	// fraction exited over all time (for steady sims)
	std::vector<double> fraction_exited_3d_surface;	// fraction exited this time step
	std::vector<double> centrelines_deposition_fraction;	// deposition fraction in each 2D/3D airways
	std::vector<double> centrelines_particle_fraction;	// deposition fraction in each 2D/3D airways

	std::vector<double> empty_vec;

	// some reading timesteps stuff
	// tells you what time step and time the read solutions are from
	unsigned int read_time_step_1;
	unsigned int read_time_step_2;
	double read_time_1;
	double read_time_2;

	bool no_motion_end_particle_sim;
	unsigned int num_wall_bdy_ids;
	double total_fraction_added;	// this is 0d only

	double particle_fraction_added_so_far;
	double particle_fraction_combined;
	double particle_fraction_0d;
	double particle_fraction_3d;
	double deposition_fraction_combined;
	double deposition_fraction_0d;
	double deposition_fraction_3d;
	double deposition_fraction_combined_max;
	double deposition_fraction_0d_max;
	double deposition_fraction_3d_max;
	double deposition_fraction_sed_0d;
	double deposition_fraction_sed_0d_max;
	double deposition_fraction_imp_0d;
	double deposition_fraction_imp_0d_max;
	double deposition_fraction_dif_0d;
	double deposition_fraction_dif_0d_max;
	double terminal_exit_fraction;
	double terminal_exit_fraction_max;

	std::vector<unsigned int> terminal_airway_to_airway_map;

	// preconditioner stuff
	PetscMatrix<Number>* pressure_mass_matrix;
	PetscMatrix<Number>* scaled_pressure_mass_matrix;
	PetscMatrix<Number>* pressure_laplacian_matrix;
	PetscMatrix<Number>* pressure_convection_diffusion_matrix;
	PetscMatrix<Number>* velocity_mass_matrix;
	PetscMatrix<Number>* velocity_matrix;

	// moghadam matrices
	/*
	PetscMatrix<Number>* moghadam_velocity_matrix;
	PetscMatrix<Number>* moghadam_velocity_preconditioner_matrix;
	PetscVector<Number>* moghadam_velocity_rhs;
	PetscVector<Number>* moghadam_velocity_temp_solution;
	*/

	SparseMatrix<Number>* moghadam_velocity_matrix;
	SparseMatrix<Number>* moghadam_velocity_preconditioner_matrix;
	SparseMatrix<Number>* moghadam_b_matrix;
	SparseMatrix<Number>* moghadam_bt_matrix;
	SparseMatrix<Number>* moghadam_c_matrix;

	SparseMatrix<Number>* moghadam_velocity_h_inverse_matrix;

	NumericVector<Number>* moghadam_velocity_rhs;
	NumericVector<Number>* moghadam_velocity_temp_solution;
	NumericVector<Number>* moghadam_schur_rhs;
	NumericVector<Number>* moghadam_schur_temp_solution;
	NumericVector<Number>* moghadam_velocity_rhs_2;
	NumericVector<Number>* moghadam_velocity_temp_solution_2;

	NumericVector<Number>* moghadam_velocity_row_sum;

	NumericVector<Number>* moghadam_lhs;
	NumericVector<Number>* moghadam_rhs;

	// moghadam schur pc
	KSP moghadam_schur_ksp;
	Mat moghadam_schur_matrix;
	MoghadamSchurShellMatrixCtx moghadam_schur_mat_ctx;

	// moghadam velocity pc
	MoghadamVelocityPCCtx  *moghadam_velocity_pc_ctx;    /* user-defined preconditioner context */

	// jacobi preconditioning
	NumericVector<Number>* diagonal;
	NumericVector<Number>* velocity_diagonal;
	NumericVector<Number>* pressure_diagonal;

	// cause it makes things easier..
	Mat schur_complement_approx;

	// Mats for submatrices of SIMPLE-type preconditioners
	IS velocity_is;
	IS pressure_is;
	Vec non_zero_cols;
	Vec non_zero_rows;
	std::vector<unsigned int> non_zero_cols_vec;
	std::vector<unsigned int> non_zero_rows_vec;
	NSShellPC  *shell;    /* user-defined preconditioner context */
	MonoShellPC  *mono_shell;    /* user-defined preconditioner context */
	NSShellPC  *schur_stokes_shell;    /* user-defined preconditioner context */
	NSShellPC  *post_solve_shell;    /* user-defined preconditioner context */
	SIMPLEShellPC  *simple_shell;    /* user-defined preconditioner context */

	PCD2ShellMatrixCtx mat_ctx;
	MonolithicMonitorCtx *mono_ctx;
	InitialGuessCtx *initial_guess_ctx;


	//ExodusII_IO_Extended exo_3d;
	//ExodusII_IO_Extended exo_1d;


};

// **************************************************************** //

