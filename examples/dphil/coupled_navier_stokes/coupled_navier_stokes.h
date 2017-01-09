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


// Needed for cast to petsc matrix
#include "libmesh/petsc_matrix.h"
#include "libmesh/petsc_vector.h"



#include "libmesh/point_locator_base.h"
#include "libmesh/point_locator_base.h"

// for element type names
#include "libmesh/elem_type.h"
#include "libmesh/string_to_enum.h"

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
extern PetscErrorCode Monolithic2ShellPCSetUp(PC,Mat,Vec,Vec,KSP,KSP);
extern PetscErrorCode MonolithicShellPCApply(PC,Vec x,Vec y);
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

	void output_command_line_options();

	void set_auto_fieldsplit_parameters();

	void read_particle_parameters();

	void output_particle_parameters();

	void read_timesteps();

	void setup_3d_mesh(EquationSystems* _es, Mesh& _mesh);

	void setup_3d_system(TransientLinearImplicitSystem * system);

	void prerefine_3d_mesh();

	void setup_1d_mesh ();

	void read_1d_mesh ();

	void generate_1d_mesh ();

	void setup_1d_system(TransientLinearImplicitSystem * system);

	void output_sim_data(bool header=false);

	void output_linear_iteration_count(bool header=false);

	void output_particle_data(bool header=false);

	void output_deposition_data(bool header=false);

	void print_flux_and_pressure();

	// write the solution at a specific time and time_step, 
	// backup=TRUE if backup is to be written
	void write_3d_solution();

	void write_particles();

	void write_1d_solution();

	void write_efficiency_solution();

	void update_times();

	void update_time_scaling();

	int solve_1d_system();

	void set_radii();

	void set_poiseuille();

	void set_efficiency();

	void calculate_1d_boundary_values();

	//returns true if wanna break iteration
	bool solve_3d_system_iteration(TransientLinearImplicitSystem * system);

	double solve_and_assemble_3d_system(TransientLinearImplicitSystem * system);

	void calculate_3d_boundary_values();

	void adaptively_refine();

	void end_timestep_admin();

	void add_1d_tree_to_mesh(std::vector<Point>& vertices, std::vector<std::vector<unsigned int> >& cell_vertices,unsigned int subdomain_id, unsigned int boundary_id);

	void create_1d_tree(std::vector<Point>& vertices, std::vector<std::vector<unsigned int> >& cell_vertices, unsigned int num_generations);

	void calculate_num_alveloar_generations_for_tree();

	// convert the 1d part of a vector from monomial to nodal
	void convert_1d_monomial_to_nodal(NumericVector<Number>& vector);

	// convert the 1d part of a vector from nodal to monomial
	void convert_1d_nodal_to_monomial(NumericVector<Number>& vector);

	void plot_error(ExodusII_IO_Extended& io);

	void write_elem_pid_1d(ExodusII_IO_Extended& io);

	void write_elem_pid_3d(ExodusII_IO_Extended& io);

	void set_elem_proc_id_3d();
	void set_elem_proc_id_1d();

	double set_double_parameter(GetPot _infile, std::string name, double default_value);
	unsigned int set_unsigned_int_parameter(GetPot _infile, std::string name, unsigned int default_value);
	int set_int_parameter(GetPot _infile, std::string name, int default_value);
	bool set_bool_parameter(GetPot _infile, std::string name, bool default_value);
	std::string set_string_parameter(GetPot _infile, std::string name, std::string default_value);
	std::string set_string_command_line_parameter(std::string name, std::string default_value);

	double scaled_time();

	void scale_3d_solution_vector(double velocity_scaling,double pressure_scaling);

	void scale_1d_solution_vector(double flux_scaling,double pressure_scaling);

	void init_dof_variable_vector(TransientLinearImplicitSystem * system, std::vector<int>& _dof_variable_type);

	void init_dof_variable_vectors();

	//void setup_results_vector(std::string results_folder, 
	//			std::vector<NumericVector<Number>* >& results_vector);		//setup vectors that stores the results from two previous simulations

	void output_poiseuille_resistance_per_generation();
	//okay want some functions that output scaled quantities
	//we can do this by constructing temporary explicit equation systems?

	void init_particles();

	void move_particles();


	void setup_variable_scalings_3D();

	void setup_variable_scalings_1D();

	void output_coupling_points();



	void construct_petsc_options_string();

	int setup_preconditioners(TransientLinearImplicitSystem * system, Mat& B);

	int compute_and_output_eigenvalues(TransientLinearImplicitSystem * system);

	int setup_is_simple (TransientLinearImplicitSystem * sys, PC my_pc);

	int construct_schur_stokes_matrix(TransientLinearImplicitSystem * sys, KSP schur_ksp);



	int test_post_solve(TransientLinearImplicitSystem * sys);

	// read the input boundary conditions for an uncoupled simulation
	int read_input_boundary_conditions();

	void calculate_1d_linear_resistance_values();

	void update_3d_dirichlet_boundary_conditions();

	void petsc_clean_up();


private:

	unsigned int sim_type;	//0 - 3d, 1 - 1d, 2 - uncoupled, 3 - expl coupled
	std::vector<std::vector<double> > element_data;
	std::vector<Airway> airway_data;
	std::vector<std::vector<double> > centreline_element_data;	//element data of centrelines
	std::vector<Airway> centreline_airway_data;
	Mesh mesh;
	MeshData mesh_data;
	MeshRefinement mesh_refinement;
	AutoPtr<EquationSystems> es;
	TransientLinearImplicitSystem * system_3d;
	TransientLinearImplicitSystem * system_1d;
	TransientLinearImplicitSystem * system_coupled;
	TransientLinearImplicitSystem * system_neumann;
	ExplicitSystem * extra_1d_data_system;
	ExplicitSystem * extra_3d_data_system;
	ExplicitSystem * particle_deposition_system_1d;
	unsigned int t_step;	
	double time;
	double dt;
	double previous_dt;
	AutoPtr<NSAssembler3D> picard;	//instead of being picard, this should be 3d assembly function
	AutoPtr<NavierStokesAssembler> ns_assembler;
	AutoPtr<CoupledAssembler> coupled_assembler;
	AutoPtr<ErrorEstimator> error_estimator;
	ErrorVector error_vector;
	std::vector<double> re_vec;
	bool reduce_dt;
	bool increase_dt;
	bool refine_mesh;
	unsigned int steps_since_last_dt_change;
	unsigned int unsteady;
	bool restart;
	PerfLog perf_log;
	AutoPtr<NumericVector<Number> > last_nonlinear_soln;	//could be pointer but too much admin
	AutoPtr<NumericVector<Number> > old_global_solution;	//could be pointer but too much admin
	Real previous_nonlinear_residual;
	Real current_nonlinear_residual;
	bool sim_3d;
	bool sim_1d;
	std::ofstream output_file;
	std::ofstream eigenvalues_file;
	std::ofstream parameter_output_file;
	std::ofstream particle_output_file;
	std::ofstream deposition_output_file;
	std::ofstream linear_iterations_output_file;
	std::vector<unsigned int> boundary_ids;
	// these are not boundary conditions but values on the respective domains
	// calculated after the fact.
	// e.g. may want to apply 3d values as boundary conditions
	std::vector<double> pressure_values_3d;
	std::vector<double> flux_values_3d;
	std::vector<double> pressure_values_1d;
	std::vector<double> flux_values_1d;
	std::vector<double> previous_pressure_values_3d;
	std::vector<double> previous_flux_values_3d;
	std::vector<double> previous_previous_flux_values_3d;
	std::vector<double> previous_pressure_values_1d;
	std::vector<double> previous_flux_values_1d;
	std::vector<double> input_pressure_values_3d;
	std::vector<double> input_flux_values_1d;
	std::vector<double> linear_resistance_values_1d;
	AutoPtr<AugmentSparsityOnInterface> augment_sparsity;
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
	//std::vector<NumericVector<Number>* > results_vector_1;		//could be a vector of pointers but too much admin
	//std::vector<NumericVector<Number>* > results_vector_2;		//could be a vector of pointers but too much admin
	bool threed;
	std::vector<unsigned int> subdomains_1d;		//list of subdomain ids that are 1d (i.e. 1d airways)
	std::vector<unsigned int> subdomains_2d;		//list of subdomain ids that are 2d (i.e. neumann boundary)
	std::vector<unsigned int> subdomains_3d;		//list of subdomain ids that are 3d (i.e. 3d airways) 
	unsigned int total_nonlinear_iterations;		//list of subdomain ids that are 3d (i.e. 3d airways)
	unsigned int local_linear_iterations;		//list of subdomain ids that are 3d (i.e. 3d airways) 
	int total_linear_iterations;
	int total_max_iterations;
	//SurfaceBoundary inflow_surface_boundary_object;
	std::vector<SurfaceBoundary* > surface_boundaries;
	unsigned int nonlinear_iteration;
	unsigned int num_linear_iterations;
	//edge num is idx, 0-elem_num 1-node1 2-node2 3-generation 4-order
	std::vector<std::vector<double> > edge_1d_data;
	unsigned int num_edges_1d;
	//node num is idx, 0-xcoord 1-ycoord 2-zcoord 3-radius
	std::vector<std::vector<double> > node_1d_data;
	unsigned int num_nodes_1d;
	std::vector<unsigned int> num_generations; //num_generations for each tree
	bool ic_set;

	// particle deposition stuff
	unsigned int particle_deposition;
	std::vector<double> timestep_sizes;	//vector of timestep sizes, 1st timestep is 0th element obv
	std::vector<Particle> particles_3D;

	std::vector<double> var_scalings_3D;
	std::vector<double> var_scalings_1D;
	std::vector<double> total_efficiency;

	std::vector<int> subtree_starting_elements;	//ending elements of the 3d mesh and the subtree that they go to
	std::vector<int> boundary_id_to_tree_id;
	std::vector<int> tree_id_to_boundary_id;
	std::vector<double> coupling_pressure_values_1d;
	std::vector<double> coupling_flux_values_1d;

	std::vector<int> centreline_terminal_id_to_tree_id;
	std::vector<Point> centreline_points_1;
	std::vector<Point> centreline_points_2;

	// this is a vector of the coupling points for each tree
	std::vector<Point> coupling_points;

	std::vector<std::vector<double> > alveolar_efficiency_per_generation;	// [timestep][generation] alveolar efficiency of deposition summed over all branches in this generation
	std::vector<std::vector<double> > tb_efficiency_per_generation;	// [timestep][generation] tracheo-bronchial efficiency of deposition summed over all branches in this generation
	int total_particles_inhaled;
	int stokes_gmres_iterations;
	bool shell_pc_created;
	bool schur_shell_pc_created;
	bool mono_shell_pc_created;

	std::vector<std::vector<double> > all_input_pressure_values_3d;
	std::vector<std::vector<double> > all_input_flux_values_0d;
	std::string input_boundary_conditions_filename;

	// preconditioner stuff
	PetscMatrix<Number>* pressure_mass_matrix;
	PetscMatrix<Number>* scaled_pressure_mass_matrix;
	PetscMatrix<Number>* pressure_laplacian_matrix;
	PetscMatrix<Number>* pressure_convection_diffusion_matrix;
	PetscMatrix<Number>* velocity_mass_matrix;
	PetscMatrix<Number>* velocity_matrix;

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


};

// **************************************************************** //

