#include "coupled_navier_stokes.h"

// ********************************************************** //
//
// This file includes all the PETSc Shell type functions
//
// setup_3d_mesh
// setup_3d_system
// prerefine_3d_mesh
// calculate_3d_boundary_values
// write_3d_solution
// solve_3d_system_iteration
// adaptively_refine
// plot_error
//
// *********************************************************** //

 

PetscErrorCode custom_outer_monitor(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy)
{
	PetscErrorCode ierr;

	//get the pc
	PC pc;
	ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);

	//get the inner iteration ksp
	int num_splits = 0;
	KSP* subksp;	//subksps
	ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);

	int num_outer_its;
	int num_inner_velocity_its;
	ierr = KSPGetIterationNumber(ksp,&num_outer_its); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(subksp[0],&num_inner_velocity_its); CHKERRQ(ierr);

	//ierr = PetscPrintf(PETSC_COMM_WORLD,"  in ksp custom monitor: num_velocity_its = %D\n",num_inner_pressure_its); CHKERRQ(ierr);

	MonolithicMonitorCtx *mono_ctx = (MonolithicMonitorCtx*)dummy;
	mono_ctx->total_velocity_iterations += num_inner_velocity_its;
	//mono_ctx->total_convection_diffusion_iterations = 1000;
	//total_inner_pressure_its += num_inner_pressure_its;
  return(0);
}


PetscErrorCode custom_inner_pressure_monitor(KSP ksp,PetscInt n,PetscReal rnorm,void *dummy)
{
	PetscErrorCode ierr;

	//get the pc
	PC pc;
	ierr = KSPGetPC(ksp,&pc); CHKERRQ(ierr);

	//get the inner iteration ksp
	int num_splits = 0;
	KSP* subksp;	//subksps
	ierr = PCFieldSplitGetSubKSP(pc,&num_splits,&subksp); CHKERRQ(ierr);

	int num_outer_its;
	int num_inner_velocity_its;
	int num_inner_pressure_its;
	ierr = KSPGetIterationNumber(ksp,&num_outer_its); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(subksp[0],&num_inner_velocity_its); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(subksp[1],&num_inner_pressure_its); CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"    in ksp custom pressure monitor: num_velocity_its = %D\n",num_inner_velocity_its); CHKERRQ(ierr);

	//total_inner_velocity_its += num_inner_velocity_its;
  return(0);
}



// ***** PCSHELL preconditioner
/*
   SampleShellPCCreate - This routine creates a user-defined
   preconditioner context.

   Output Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode ShellPCCreate(NSShellPC **shell)
{
  NSShellPC  *newctx;
  PetscErrorCode ierr;

  //ierr         = PetscNew(NSShellPC,&newctx);CHKERRQ(ierr);
	ierr         = PetscNew(&newctx);CHKERRQ(ierr);
  newctx->pressure_mass_matrix = 0;
  newctx->pressure_laplacian_matrix = 0;
	newctx->pressure_laplacian_preconditioner = 0;
  newctx->pressure_convection_diffusion_matrix = 0;
  newctx->velocity_matrix = 0;
  newctx->S_approx = 0;
  newctx->inner_mass_ksp = 0;
  newctx->inner_lap_ksp = 0;
  newctx->inner_velocity_ksp = 0;
  newctx->temp_vec = 0;
  newctx->temp_vec_2 = 0;
  *shell       = newctx;
  return 0;
}


// ***** PCSHELL preconditioner
/*
   SampleShellPCCreate - This routine creates a user-defined
   preconditioner context.

   Output Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode ShellPCCreate(SIMPLEShellPC **shell)
{
  SIMPLEShellPC  *newctx;
  PetscErrorCode ierr;

  //ierr         = PetscNew(NSShellPC,&newctx);CHKERRQ(ierr);
	ierr         = PetscNew(&newctx);CHKERRQ(ierr);
  *shell       = newctx;
  return 0;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleShellPCSetUp"
/*
   SampleShellPCSetUp - This routine sets up a user-defined
   preconditioner context.

   Input Parameters:
.  pc    - preconditioner object
.  pmat  - preconditioner matrix
.  x     - vector

   Output Parameter:
.  shell - fully set up user-defined preconditioner context

   Notes:
   In this example, we define the shell preconditioner to be Jacobi's
   method.  Thus, here we create a work vector for storing the reciprocal
   of the diagonal of the preconditioner matrix; this vector is then
   used within the routine SampleShellPCApply().
*/
PetscErrorCode PressureShellPCSetUp(PC pc,Mat pressure_mass_matrix, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);
  //Vec            diag;

	//ierr = KSPDestroy(&shell->inner_lap_ksp);CHKERRQ(ierr);
	//ierr = VecDestroy(&shell->temp_vec);CHKERRQ(ierr);
	//ierr = VecDestroy(&shell->temp_vec_2);CHKERRQ(ierr);


  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SET MASS MATRIX ******************* //

  shell->pressure_mass_matrix = pressure_mass_matrix;

	// ********* SETUP MASS MATRIX KSP **************** //

	//ierr = KSPDestroy(&shell->inner_mass_ksp);CHKERRQ(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_mass_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_mass_ksp,Amat,Pmat); CHKERRQ(ierr);

	// setup up the operators for the first inner solve

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);

	// set the pc to be a ksp
	// should probably use gmres by default okay
	ierr = PCSetType(local_mass_pc,PCKSP); CHKERRQ(ierr);
	KSP local_mass_ksp;

	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(local_mass_ksp,"mass_matrix_"); CHKERRQ(ierr);
	// setup up the operators for the first inner solve
	ierr = KSPSetOperators(local_mass_ksp,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_mass_pc,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_mass_ksp); CHKERRQ(ierr);

	

	ierr = PetscLogStagePop();

//	printf ("inside shell pc setup");
  return 0;
}

PetscErrorCode PCDShellPCSetUp(PC pc,Mat pressure_mass_matrix,Mat pressure_laplacian_matrix,Mat pressure_laplacian_preconditioner,Mat pressure_convection_diffusion_matrix, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP MATRICES ****************** //

	shell->pressure_mass_matrix = pressure_mass_matrix;
  shell->pressure_laplacian_matrix = pressure_laplacian_matrix;
	shell->pressure_laplacian_preconditioner = pressure_laplacian_preconditioner;
  shell->pressure_convection_diffusion_matrix = pressure_convection_diffusion_matrix;


	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_lap_ksp,Amat,Pmat); CHKERRQ(ierr);


	PC local_lap_pc;
	ierr = KSPGetPC(shell->inner_lap_ksp,&local_lap_pc); CHKERRQ(ierr);

	// set the pc to be a ksp
	// should probably use gmres by default okay
	ierr = PCSetType(local_lap_pc,PCKSP); CHKERRQ(ierr);


	KSP local_lap_ksp;
	ierr = PCKSPGetKSP(local_lap_pc,&local_lap_ksp); CHKERRQ(ierr);


	ierr = KSPSetOptionsPrefix(local_lap_ksp,"laplacian_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_lap_ksp); CHKERRQ(ierr);
	

	// setup up the operators for the first inner solve
	ierr = KSPSetOperators(local_lap_ksp,shell->pressure_laplacian_matrix,shell->pressure_laplacian_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_lap_pc,shell->pressure_laplacian_matrix,shell->pressure_laplacian_matrix); CHKERRQ(ierr);

	


	// ********* SETUP MASS KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_mass_ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_mass_ksp,Amat,Pmat); CHKERRQ(ierr);

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);
	ierr = PCSetType(local_mass_pc,PCKSP); CHKERRQ(ierr);
	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);

	
	ierr = KSPSetOptionsPrefix(local_mass_ksp,"mass_matrix_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_mass_ksp); CHKERRQ(ierr);
	
	ierr = KSPSetOperators(local_mass_ksp,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_mass_pc,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);




  // ******************** SETUP VECTORS ***************************** //
	PetscInt n_col;
	PetscInt n_col_local;
	ierr = MatGetSize(Amat,NULL,&n_col); CHKERRQ(ierr);
	ierr = MatGetLocalSize(Amat,NULL,&n_col_local); CHKERRQ(ierr);
	// create some vectors that can hold temporary info
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec_2); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec_2, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec_2); CHKERRQ(ierr);


//	printf ("inside shell pc setup");

	ierr = PetscLogStagePop();
  return 0;
}


PetscErrorCode PCD2ShellPCSetUp(PC pc,Mat pressure_mass_matrix,Mat pressure_laplacian_matrix,Mat pressure_laplacian_preconditioner,Mat pressure_convection_diffusion_matrix, KSP schur_ksp)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);

  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	//create a KSP that can do the solving we require
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_mass_ksp); CHKERRQ(ierr);
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	//MatStructure mat_structure;
	//ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat,&mat_structure); CHKERRQ(ierr);
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOperators(shell->inner_mass_ksp,Amat,Pmat); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_lap_ksp,Amat,Pmat); CHKERRQ(ierr);


  // create some vectors for temporary pc applying
	PetscInt n_col;
	PetscInt n_col_local;
	ierr = MatGetSize(Amat,NULL,&n_col); CHKERRQ(ierr);
	ierr = MatGetLocalSize(Amat,NULL,&n_col_local); CHKERRQ(ierr);
	// create some vectors that can hold temporary info
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&shell->temp_vec_2); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetSizes(shell->temp_vec_2, n_col_local, n_col); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec); CHKERRQ(ierr);
	ierr = VecSetFromOptions(shell->temp_vec_2); CHKERRQ(ierr);


	shell->pressure_mass_matrix = pressure_mass_matrix;
  shell->pressure_laplacian_matrix = pressure_laplacian_matrix;
	shell->pressure_laplacian_preconditioner = pressure_laplacian_preconditioner;
  shell->pressure_convection_diffusion_matrix = pressure_convection_diffusion_matrix;

//	printf ("inside shell pc setup");
	ierr = PetscLogStagePop();
  return 0;
}


/*
   MonolithicShellPCSetUp - This routine sets up the monolithic shell preconditioner
	when we don't construct the schur complement columnwise. mostly defunct.
*/
PetscErrorCode MonolithicShellPCSetUp(PC pc,Mat velocity_matrix, KSP schur_ksp)
{
  	NSShellPC  *shell;
  	PetscErrorCode ierr;

	//ierr = PetscLogStagePush(2); // not sure why this doesn't work anymore


  	ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SET THE MATRIX ******************* //

  	shell->velocity_matrix = velocity_matrix;

	// ********* SETUP VELOCITY MATRIX KSP **************** //

	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->velocity_matrix,shell->velocity_matrix); CHKERRQ(ierr);

	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	// set the pc to be used, lu for now
	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"ns3d1d_fieldsplit_1_pc_velocity_"); CHKERRQ(ierr);
	// setup up the operators for the first inner solve
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);



	// Some of the Mats we will extract and create
	Mat Bt;
	Mat B;
	Mat A11;
	Mat S;


	Mat Bt_dense;	// dense version of BT
	Mat T;
	Mat F;	// factored velocity matrix
	IS perm, iperm;
	MatFactorInfo info;
	Mat S_approx_dense;	// approximate dense schur complement

	ierr = KSPGetOperators(schur_ksp,&S,NULL); CHKERRQ(ierr);
	ierr = MatSchurComplementGetSubMatrices(S,NULL,NULL,&Bt,&B,&A11);
	printf ("hi\n");


	ierr = MatConvert(Bt,MATDENSE,MAT_INITIAL_MATRIX,&Bt_dense);
	ierr = MatDuplicate(Bt_dense,MAT_COPY_VALUES,&T);

	// create the factorisation
	printf ("2\n");
	
  ierr = MatGetOrdering(shell->velocity_matrix,  MATORDERINGNATURAL,  &perm,  &iperm); CHKERRQ(ierr);     
  ierr = MatFactorInfoInitialize(&info); CHKERRQ(ierr);
	
	// superlu doesn't work with matmatsolve, superlu_dist does.
	ierr = MatGetFactor(shell->velocity_matrix,MATSOLVERSUPERLU_DIST,MAT_FACTOR_LU,&F); CHKERRQ(ierr); 
	printf ("3\n");
	ierr = MatLUFactorSymbolic(F,shell->velocity_matrix,perm,iperm,&info); CHKERRQ(ierr); 
	printf ("4\n");
	ierr = MatLUFactorNumeric(F,shell->velocity_matrix,&info); CHKERRQ(ierr); 
	

	printf ("5\n");
	// if using matlufactor then use shell->velocity_matrix, is using matgetfactor use F
	ierr = MatMatSolve(F,Bt_dense,T);

	printf ("6\n");
	ierr = MatMatMult(B,T,MAT_INITIAL_MATRIX,1.0,&S_approx_dense);
	printf ("7\n");

	// convert S_approx to sparse
	ierr = MatConvert(S_approx_dense,MATAIJ,MAT_INITIAL_MATRIX,&shell->S_approx);

	// S_approx = A11 - S_approx
	ierr = MatAXPY(shell->S_approx,-1.0,A11,DIFFERENT_NONZERO_PATTERN);
	ierr = MatScale(shell->S_approx,-1.0);


	// only use A11, then we don't need the matscale below
	//ierr = MatDuplicate(A11,MAT_COPY_VALUES,&shell->S_approx);
	printf ("8\n");

	// setup up the operators for the 
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);

	printf ("yeah\n");

	// delete the matrices we have created
	ierr = MatDestroy(&Bt_dense); CHKERRQ(ierr);
	ierr = MatDestroy(&T); CHKERRQ(ierr);
	
	ierr = MatDestroy(&F); CHKERRQ(ierr);
	ierr = ISDestroy(&perm); CHKERRQ(ierr);
	ierr = ISDestroy(&iperm); CHKERRQ(ierr);
	ierr = MatDestroy(&S_approx_dense); CHKERRQ(ierr);
	

	//ierr = PetscLogStagePop();

	printf ("inside shell pc monolithic velocity setup");
  return 0;
}


/*
   Monolithic2ShellPCSetUp - This solves the matrix in the schur complement vector by vector.
*/
PetscErrorCode Monolithic2ShellPCSetUp(PC pc,Mat velocity_matrix, Vec non_zero_cols, Vec non_zero_rows, KSP schur_ksp)
{
  	NSShellPC  *shell;
  	PetscErrorCode ierr;

	//ierr = PetscLogStagePush(2);	// not sure why this log thing makes an error now...

  	ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SET MASS MATRIX ******************* //

  	shell->velocity_matrix = velocity_matrix;

	// ********* SETUP VELOCITY MATRIX KSP **************** //

	// create
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);
	
	// set the operators for the velocity matrix inversion
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->velocity_matrix,shell->velocity_matrix); CHKERRQ(ierr);

	// setup up the PC for the velocity matrix inversion
	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	// set the pc to be used, lu for now
	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"ns3d1d_fieldsplit_1_pc_velocity_"); CHKERRQ(ierr);
	// setup up the operators for the first inner solve
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);


	// ************* Figure out how many 


	// Get the matrices we need from the system
	Mat Bt;
	Mat B;
	Mat A11;
	Mat S;

	ierr = KSPGetOperators(schur_ksp,&S,NULL); CHKERRQ(ierr);
	ierr = MatSchurComplementGetSubMatrices(S,NULL,NULL,&Bt,&B,&A11);

	// Create the Factor matrix object and associated items
	Mat F;	// factored velocity matrix
	IS perm, iperm;
	MatFactorInfo info;

	printf ("hi\n");


	printf ("2\n");
	
	// initialise the items for the factorisation (probably unnecessary)
  	ierr = MatGetOrdering(shell->velocity_matrix,  MATORDERINGNATURAL,  &perm,  &iperm); CHKERRQ(ierr);     
  	ierr = MatFactorInfoInitialize(&info); CHKERRQ(ierr);

	// Do the factorisation
	// superlu doesn't work with matmatsolve, superlu_dist does.
	// can use different arguments in here for the solver
	// 1) MATSOLVERPETSC - very slow
	// 2) MATSOLVERSUPERLU_DIST - works nicely
	ierr = MatGetFactor(shell->velocity_matrix,MATSOLVERSUPERLU_DIST,MAT_FACTOR_LU,&F); CHKERRQ(ierr); 
	printf ("3\n");
	ierr = MatLUFactorSymbolic(F,shell->velocity_matrix,perm,iperm,&info); CHKERRQ(ierr); 
	printf ("4\n");
	ierr = MatLUFactorNumeric(F,shell->velocity_matrix,&info); CHKERRQ(ierr); 
	

	printf ("5\n");

	// **** MANUALLY MAKE MATRIX BF^-1BT

	// instead of matmatsolve we want to do each column and then put it into a new sparse matrix T.
	PetscInt num_coupling_points;
	ierr = VecGetSize(non_zero_cols,&num_coupling_points);
	
	PetscInt m;
	PetscInt m_local;
	PetscInt n_1d;
	PetscInt n_1d_local;
	ierr = MatGetSize(Bt,&m,&n_1d);
	ierr = MatGetLocalSize(Bt,&m_local,&n_1d_local);
	Vec Bt_col;
	VecCreate(PETSC_COMM_WORLD,&Bt_col);
	VecSetSizes(Bt_col,m_local,m);
	VecSetFromOptions(Bt_col);
	Vec T_col;
	VecCreate(PETSC_COMM_WORLD,&T_col);
	VecSetSizes(T_col,m_local,m);
	VecSetFromOptions(T_col);

	PetscInt n;
	PetscInt n_local;
	ierr = MatGetSize(B,NULL,&n);
	ierr = MatGetLocalSize(B,NULL,&n_local);

	Vec B_row;
	VecCreate(PETSC_COMM_WORLD,&B_row);
	VecSetSizes(B_row,n_local,n);
	VecSetFromOptions(B_row);

	std::cout << "n = " << n << std::endl;
	std::cout << "m = " << m << std::endl;

	// Create a matrix, T_vals of the values that will go into the BF^-1BT matrix
	// T_rows are the what rows they go into (the flux coupling dofs)
	// T_cols are the what rows they go into (the pressure coupling dofs)
	PetscScalar T_vals[num_coupling_points][num_coupling_points];
	PetscInt T_rows[num_coupling_points];
	PetscInt T_cols[num_coupling_points];

	// Create a matrix to put the result into
	// We know that there won't be more than num_coupling_points on each row,
	// so use this as an estimate.
	Mat new_T;
	MatCreateAIJ(PETSC_COMM_WORLD,n_1d_local,n_1d_local,n_1d,n_1d,
                       num_coupling_points,NULL,num_coupling_points,NULL,&new_T);

	// loop over the columns
	for(unsigned int i=0; i<num_coupling_points; i++)
	{
		PetscScalar col_number[1];
		PetscInt coupling_point[1] = {i};	// the coupling point we are interested in (must be an array...)
		ierr = VecGetValues(non_zero_cols,1,coupling_point,col_number);

		// extract the correct column as a Vec
		ierr = MatGetColumnVector(Bt,Bt_col,col_number[0]);

		// copmute F^-1BT for this particular row
		ierr = MatSolve(F,Bt_col,T_col);

		// loop over the rows for the other bit of the matrix multiplication		
		for(unsigned int j=0; j<num_coupling_points; j++)
		{
			PetscScalar row_number[1];
			coupling_point[0] = j;			// the coupling point we are interested in (must be an array...)
			ierr = VecGetValues(non_zero_rows,1,coupling_point,row_number);		// get the actual dof index in the B matrix
		
			// create some arrays for the values and indices in the rows
			const PetscScalar *B_row_array;
			const PetscInt *B_row_cols;
			PetscInt B_row_n_cols;

			// extract the correct row
			ierr = MatGetRow(B,row_number[0],&B_row_n_cols,&B_row_cols,&B_row_array);
			VecSet(B_row,0);	// reinitialise each time

			// convert the row values to a Vec for dot product
			ierr = VecSetValues(B_row,B_row_n_cols,B_row_cols,B_row_array,INSERT_VALUES);

			// do a dot product with the column that was create in the outer loop
			PetscScalar T_val;
			ierr = VecDot(B_row,T_col,&T_val);

			// save the result in the correct place and save the row and column indices
			T_vals[i][j] = T_val;
			T_cols[i] = col_number[0];
			T_rows[j] = row_number[0];

			
		}
		
	}

	// set the values of the matrix and assemble them
	MatSetValues(new_T,num_coupling_points,T_rows,num_coupling_points,T_cols,(const PetscScalar*)T_vals,INSERT_VALUES);
	MatAssemblyBegin(new_T,MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd(new_T,MAT_FINAL_ASSEMBLY );

	// put the values into the correct matrix... for no good reason..
	ierr = MatDuplicate(new_T,MAT_COPY_VALUES,&shell->S_approx);
	
	// S_approx = A11 - S_approx
	ierr = MatAXPY(shell->S_approx,-1.0,A11,DIFFERENT_NONZERO_PATTERN);
	ierr = MatScale(shell->S_approx,-1.0);

	printf ("8\n");

	// setup up the operators for the solve later
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);

	printf ("yeah\n");

	// delete the matrices we have created
	ierr = VecDestroy(&Bt_col);
	ierr = VecDestroy(&T_col);
	ierr = VecDestroy(&B_row);
	ierr = MatDestroy(&new_T);
	ierr = MatDestroy(&F); CHKERRQ(ierr);
	ierr = ISDestroy(&perm); CHKERRQ(ierr);
	ierr = ISDestroy(&iperm); CHKERRQ(ierr);
	

	//ierr = PetscLogStagePop();

	printf ("inside shell pc monolithic 2 velocity setup");
  return 0;
}




/*
   Monolithic3ShellPCSetUp - This solves the matrix in the schur complement vector by vector.
*/
PetscErrorCode Monolithic3ShellPCSetUp(PC pc,Mat schur_complement_approx, KSP schur_ksp)
{
  	NSShellPC  *shell;
  	PetscErrorCode ierr;

	//ierr = PetscLogStagePush(2);	// not sure why this log thing makes an error now...

  	ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP VELOCITY MATRIX KSP **************** //

	// create
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);

	// set the pc to be used, lu for now
	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"ns3d1d_fieldsplit_1_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);




	// Get the matrices we need from the system
	Mat A11;
	Mat S;

	ierr = KSPGetOperators(schur_ksp,&S,NULL); CHKERRQ(ierr);
	ierr = MatSchurComplementGetSubMatrices(S,NULL,NULL,NULL,NULL,&A11);

	// copy the schur_complement_approx containing the stokes approximation to the S_approx
	ierr = MatDuplicate(schur_complement_approx,MAT_COPY_VALUES,&shell->S_approx);
	
	// S_approx = A11 - S_approx
	ierr = MatAXPY(shell->S_approx,-1.0,A11,DIFFERENT_NONZERO_PATTERN);
	ierr = MatScale(shell->S_approx,-1.0);






	// setup up the operators for the solve later
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);	

	//ierr = PetscLogStagePop();

	printf ("inside shell pc monolithic 3 velocity setup");
  return 0;
}






PetscErrorCode LSCShellPCSetUp(PC pc, KSP schur_ksp)
{
	printf ("inside shell pc lsc james setup");
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_lap_ksp,"laplacian_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_lap_ksp); CHKERRQ(ierr);
	



	
	// ********** SETUP THE INTERMEDIATE VECS ***************	//
	Mat A,Ap,B,C;
	MatSchurComplementGetSubMatrices(Amat,&A,&Ap,&B,&C,NULL);
	MatGetVecs(A,&shell->temp_vec,&shell->temp_vec_2);
	MatGetVecs(Amat,&shell->temp_vec_3,NULL);
	VecDuplicate(shell->temp_vec,&shell->lsc_scale);

	// ********** SETUP L = A10 * A01 *************************	//
	MatMatMult(C,B,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&shell->lsc_laplacian_matrix);		// could possibly use MAT_REUSE_MATRIX

	// ********** SETUP DIAGONAL ******************* //
	MatGetDiagonal(Ap,shell->lsc_scale); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(shell->lsc_scale);
	
	// setup up the operators for the laplacian solve
	ierr = KSPSetOperators(shell->inner_lap_ksp,shell->lsc_laplacian_matrix,shell->lsc_laplacian_matrix); CHKERRQ(ierr);


	printf ("done shell pc lsc james setup");

	ierr = PetscLogStagePop();
  return 0;
}







PetscErrorCode LSCScaledShellPCSetUp(PC pc, Mat velocity_mass_matrix, KSP schur_ksp)
{
	printf ("inside shell pc lsc james setup");
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_lap_ksp,"laplacian_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_lap_ksp); CHKERRQ(ierr);
	

	// ********* SET THE VELOCITY MASS MATRIX ****************	//
  	shell->velocity_matrix = velocity_mass_matrix;


	
	// ********** SETUP THE INTERMEDIATE VECS ***************	//
	Mat A,Ap,B,C,B_scaled;
	MatSchurComplementGetSubMatrices(Amat,&A,&Ap,&B,&C,NULL);
	MatGetVecs(A,&shell->temp_vec,&shell->temp_vec_2);
	MatGetVecs(Amat,&shell->temp_vec_3,NULL);
	VecDuplicate(shell->temp_vec,&shell->lsc_scale);
	MatDuplicate(B,MAT_COPY_VALUES,&B_scaled);	// create temporary place for B_scaled

	// ********** SETUP DIAGONAL ******************* //
	MatGetDiagonal(shell->velocity_matrix,shell->lsc_scale); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(shell->lsc_scale);
	
	// ********** SETUP L = A10* Q_V^-1* A01 *************************	//
	MatDiagonalScale(B_scaled,shell->lsc_scale,NULL);
	MatMatMult(C,B_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&shell->lsc_laplacian_matrix);		// could possibly use MAT_REUSE_MATRIX

	// setup up the operators for the laplacian solve
	ierr = KSPSetOperators(shell->inner_lap_ksp,shell->lsc_laplacian_matrix,shell->lsc_laplacian_matrix); CHKERRQ(ierr);

	//std::cout << "Bt" << std::endl;
	//MatView(B,PETSC_VIEWER_STDOUT_SELF);

	//std::cout << "B" << std::endl;
	//MatView(C,PETSC_VIEWER_STDOUT_SELF);


	//std::cout << "lsc_laplacian after" << std::endl;
	//MatView(shell->lsc_laplacian_matrix,PETSC_VIEWER_STDOUT_SELF);

	MatDestroy(&B_scaled);
	printf ("done shell pc lsc james setup");

	ierr = PetscLogStagePop();
  return 0;
}






PetscErrorCode LSCScaledStabilisedShellPCSetUp(PC pc, Mat velocity_mass_matrix, KSP schur_ksp)
{
	printf ("Inside shell pc lsc scaled stabilised james setup\n");
  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);



	std::cout << "hello?" << std::endl;

  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);

	// ********* SETUP LAPLACIAN KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_lap_ksp); CHKERRQ(ierr);

	//get operators from the outer pcshell solver
	Mat Amat, Pmat;
	ierr = KSPGetOperators(schur_ksp,&Amat,&Pmat); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_lap_ksp,"laplacian_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_lap_ksp); CHKERRQ(ierr);
	

	// ********* SET THE VELOCITY MASS MATRIX ****************	//
  	shell->velocity_matrix = velocity_mass_matrix;

	std::cout << "hello?" << std::endl;

	
	// ********** SETUP THE INTERMEDIATE VECS ***************	//
	Mat A,Ap,Bt,B,Bt_scaled,C;
	MatSchurComplementGetSubMatrices(Amat,&A,&Ap,&Bt,&B,&C);
	MatGetVecs(A,&shell->temp_vec,&shell->temp_vec_2);
	MatGetVecs(Amat,&shell->temp_vec_3,NULL);
	VecDuplicate(shell->temp_vec,&shell->lsc_scale);
	MatDuplicate(Bt,MAT_COPY_VALUES,&Bt_scaled);	// create temporary place for B_scaled

	std::cout << "hello?" << std::endl;
	// ********** SETUP VELOCITY MASS DIAGONAL INVERSE ******************* //
	MatGetDiagonal(shell->velocity_matrix,shell->lsc_scale); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(shell->lsc_scale);
	
	std::cout << "hello?" << std::endl;
	// ********** CALCULATE D_inv = inv(diag(B*diag(F)^-1*B^T - C)) ************ //
	Mat D;	// D matrix
	Mat D_temp;	// D matrix
	Mat B_F_inv_Bt_D_inv;	// D matrix
	Vec F_diag_inv;
	MatDuplicate(Bt,MAT_COPY_VALUES,&D_temp);	// D = B^T
	MatGetVecs(A,&F_diag_inv,NULL);
	MatGetDiagonal(A,F_diag_inv); 			// get the f_diag
	VecReciprocal(F_diag_inv);			// invert it
	MatDiagonalScale(D_temp,F_diag_inv,NULL);	// D = diag(F)^-1 * B^T
	MatMatMult(B,D_temp,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&D);	// D = B * diag(F)^-1 * B^T
	MatDuplicate(D,MAT_COPY_VALUES,&B_F_inv_Bt_D_inv);	// create temporary place for B_scaled
	MatAXPY(D,-1.0,C,DIFFERENT_NONZERO_PATTERN);	// D = B * diag(F)^-1 * B^T - C
	MatGetVecs(D,&shell->lsc_stab_alpha_D_inv,NULL);
	MatGetDiagonal(D,shell->lsc_stab_alpha_D_inv); // D = diag(B * diag(F)^-1 * B^T - C)
	VecReciprocal(shell->lsc_stab_alpha_D_inv);	// invert it

	// make B_F_inv_Bt_D_inv
	MatDiagonalScale(B_F_inv_Bt_D_inv,NULL,shell->lsc_stab_alpha_D_inv);	// B_F_inv_Bt_D_inv (the shell var doesn't have alpha yet...)
	
	
	//std::cout << "hello?" << std::endl;

	// ********** SETUP L = A10* Q_V^-1* A01 *************************	//
	Vec B_Q_inv_Bt_diag;
	Vec C_diag_inv;
	Vec D_r_vec_sqrt;
	MatDiagonalScale(Bt_scaled,shell->lsc_scale,NULL);
	MatMatMult(B,Bt_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&shell->lsc_laplacian_matrix);		// could possibly use MAT_REUSE_MATRIX

	//std::cout << "lsc_laplacian" << std::endl;
	//MatView(shell->lsc_laplacian_matrix,PETSC_VIEWER_STDOUT_SELF);
	
	
	
	//std::cout << "hello?" << std::endl;

	// ********* CALCULATE D_r_vec ***************** //
	PetscReal D_r_infty;
	//std::cout << "hi?" << std::endl;
	MatGetVecs(shell->lsc_laplacian_matrix,&B_Q_inv_Bt_diag,NULL);
	MatGetDiagonal(shell->lsc_laplacian_matrix,B_Q_inv_Bt_diag); // 
	MatGetVecs(C,&C_diag_inv,NULL);
	MatGetVecs(C,&D_r_vec_sqrt,NULL);
	MatGetDiagonal(C,C_diag_inv); // 
	VecReciprocal(C_diag_inv);	// invert it
	VecPointwiseMult(D_r_vec_sqrt,B_Q_inv_Bt_diag,C_diag_inv);
	VecNorm(D_r_vec_sqrt,NORM_INFINITY,&D_r_infty);
	VecSqrtAbs(D_r_vec_sqrt);


	//std::cout << "hello?" << std::endl;
	// we need to calculate the spectral radius, largest eigenvalue of Q_inv * F
	Mat Q_inv_F;
	MatDuplicate(A,MAT_COPY_VALUES,&Q_inv_F);	// create temporary place for B_scaled
	MatDiagonalScale(Q_inv_F,shell->lsc_scale,NULL);	// D = diag(F)^-1 * B^T



	
	// ************ do power series to calculate spectral radius ************ //
	Vec b_k;
	Vec b_k_plus_1;
	// initialise b_k
	MatGetVecs(A,&b_k,NULL);
	MatGetVecs(A,&b_k_plus_1,NULL);
	//std::cout << "hi?" << std::endl;
	MatGetDiagonal(A,b_k); // get the f_diag  
	PetscRandom rctx;   
	PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
	PetscRandomSetFromOptions(rctx);
     	VecSetRandom(b_k,rctx);	
	//std::cout << "hi?" << std::endl;
	//std::cout << "hi?" << std::endl;
	PetscReal norm = 1.0;
	PetscReal max_eigenvalue_Q_inv_F = 1.0;
	PetscReal max_eigenvalue_Q_inv_F_prev = 1.0;
	
	//std::cout << "hello?" << std::endl;
	
	for(int i=0; i<100; i++)
	{
		MatMult(Q_inv_F,b_k,b_k_plus_1);	//Ab_k = b_k+1		

		max_eigenvalue_Q_inv_F_prev = max_eigenvalue_Q_inv_F;
		// eigenvalue estimate
		VecDot(b_k,b_k_plus_1,&max_eigenvalue_Q_inv_F); // mu = b_k*A*b_k/b_k*b_k
		//max_eigenvalue = max_eigenvalue/norm;
		//std::cout << "1/norm " << i << ": " << norm << std::endl;

		// calc new norm
		VecNorm(b_k_plus_1,NORM_2,&norm);
		VecScale(b_k_plus_1,1./norm);
		VecCopy(b_k_plus_1,b_k);
	}

	std::cout << "max eigenvalue Q inv F : " << max_eigenvalue_Q_inv_F << std::endl;
	std::cout << "error : " << fabs(max_eigenvalue_Q_inv_F - max_eigenvalue_Q_inv_F_prev)/max_eigenvalue_Q_inv_F << std::endl;


	//std::cout << "hello?" << std::endl;
	

 	VecDestroy(&b_k);
	VecDestroy(&b_k_plus_1);

	MatGetVecs(B_F_inv_Bt_D_inv,&b_k,NULL);
	MatGetVecs(B_F_inv_Bt_D_inv,&b_k_plus_1,NULL);
    
	VecSetRandom(b_k,rctx);	

	PetscReal max_eigenvalue_B_F_inv_Bt_D_inv = 1.0;
	PetscReal max_eigenvalue_B_F_inv_Bt_D_inv_prev = 1.0;
	norm = 1.0;

	for(int i=0; i<100; i++)
	{
		MatMult(B_F_inv_Bt_D_inv,b_k,b_k_plus_1);	//Ab_k = b_k+1		

		max_eigenvalue_B_F_inv_Bt_D_inv_prev = max_eigenvalue_B_F_inv_Bt_D_inv;
		// eigenvalue estimate
		VecDot(b_k,b_k_plus_1,&max_eigenvalue_B_F_inv_Bt_D_inv); // mu = b_k*A*b_k/b_k*b_k
		//max_eigenvalue = max_eigenvalue/norm;
		//std::cout << "1/norm " << i << ": " << norm << std::endl;

		// calc new norm
		VecNorm(b_k_plus_1,NORM_2,&norm);
		VecScale(b_k_plus_1,1./norm);
		VecCopy(b_k_plus_1,b_k);
	}

	std::cout << "max eigenvalue B_F_inv_Bt_D_inv : " << max_eigenvalue_B_F_inv_Bt_D_inv << std::endl;
	std::cout << "error : " << fabs(max_eigenvalue_B_F_inv_Bt_D_inv - max_eigenvalue_B_F_inv_Bt_D_inv_prev)/max_eigenvalue_B_F_inv_Bt_D_inv << std::endl;


	
	//std::cout << "hello?" << std::endl;	
 	VecDestroy(&b_k);
	VecDestroy(&b_k_plus_1);


	//std::cout << "hello?" << std::endl;
	// ************ Set gamma and alpha ************** //
	PetscReal gamma = max_eigenvalue_Q_inv_F/3.0;
	//PetscReal gamma = 1.0e+2/3.0;
	gamma = gamma/D_r_infty;
	PetscReal alpha = 1.0/max_eigenvalue_B_F_inv_Bt_D_inv;


   	PetscRandomDestroy(&rctx);
	

	//std::cout << "hello?" << std::endl;

	// ********** SETUP L = A10* Q_V^-1* A01 - gamma*D_r_sqrt* C*D_r_sqrt *************************	//
	Mat D_r_sqrt_C_D_r_sqrt;
	MatDuplicate(C,MAT_COPY_VALUES,&D_r_sqrt_C_D_r_sqrt);	// create temporary place for B_scaled
	MatDiagonalScale(D_r_sqrt_C_D_r_sqrt,D_r_vec_sqrt,D_r_vec_sqrt);	// D = diag(F)^-1 * B^T
	MatAXPY(shell->lsc_laplacian_matrix,-gamma,D_r_sqrt_C_D_r_sqrt,DIFFERENT_NONZERO_PATTERN);		// could possibly use MAT_REUSE_MATRIX
	
	//std::cout << "hello?" << std::endl;
	// ********** SETUP \alpha D_inv
	VecScale(shell->lsc_stab_alpha_D_inv,alpha);
	


	

	//std::cout << "hello?" << std::endl;
	// setup up the operators for the laplacian solve
	ierr = KSPSetOperators(shell->inner_lap_ksp,shell->lsc_laplacian_matrix,shell->lsc_laplacian_matrix); CHKERRQ(ierr);

	//std::cout << "lsc_laplacian after" << std::endl;
	//MatView(shell->lsc_laplacian_matrix,PETSC_VIEWER_STDOUT_SELF);

	// *********** DESTROY ALL THE TEMPORARY THINGS VECs and MATs ***********	//
	MatDestroy(&D);	// D matrix
	MatDestroy(&D_temp);	// D matrix
	MatDestroy(&B_F_inv_Bt_D_inv);	// D matrix
	VecDestroy(&F_diag_inv);
	VecDestroy(&B_Q_inv_Bt_diag);
	VecDestroy(&C_diag_inv);
	VecDestroy(&D_r_vec_sqrt);
	MatDestroy(&Q_inv_F);
	MatDestroy(&D_r_sqrt_C_D_r_sqrt);
	MatDestroy(&Bt_scaled);

	printf ("done shell pc lsc james setup");

	//std::cout << "hello?" << std::endl;
	ierr = PetscLogStagePop();
  return 0;
}







PetscErrorCode SIMPLEShellPCSetUp(PC pc, IS velocity_is, IS pressure_is, KSP schur_ksp)
{
	printf ("inside shell pc SIMPLE james setup");
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);


	// ************* WHAT DO WE NEED ****************** //
	// we need a ksp for the velocity solve
	// we need a ksp for the schur solve


	// ************** CREATE THE SUB MATRICES WE WILL USE ********* //
	shell->velocity_is = velocity_is;
	shell->pressure_is = pressure_is;

	Mat A;
	ierr = PCGetOperators(pc,&A,NULL); CHKERRQ(ierr);

  	MatGetSubMatrix(A,shell->velocity_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->A);
  	MatGetSubMatrix(A,shell->velocity_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->Bt);
  	MatGetSubMatrix(A,shell->pressure_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->B);
  	MatGetSubMatrix(A,shell->pressure_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->C);


	// ********* SETUP VELOCITY KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"convection_diffusion_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);

	// setup up the operators for the velocity solve
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->A,shell->A); CHKERRQ(ierr);
	
	ierr = KSPSetType(shell->inner_velocity_ksp, KSPPREONLY); CHKERRQ(ierr);

	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_velocity_pc,MATSOLVERSUPERLU_DIST);

	// ************ SETUP SCHUR MATRIX ******************** //
	
	// ********** SETUP VELOCITY  DIAGONAL INVERSE ******************* //
	Vec F_diag_inv;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(F_diag_inv);


	// ********** SETUP S_approx = C -  B* diag(F)^-1* B^T *************************	//
	Mat Bt_scaled;
	Mat B_diag_F_inv_Bt;
	MatDuplicate(shell->Bt,MAT_COPY_VALUES,&Bt_scaled);	// create temporary place for Bt_scaled
	MatDuplicate(shell->C,MAT_COPY_VALUES,&B_diag_F_inv_Bt);	// create sapce for B_diag_F_inv_Bt
	MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);	// set values of S_approx
	MatDiagonalScale(Bt_scaled,F_diag_inv,NULL);

	MatMatMult(shell->B,Bt_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&B_diag_F_inv_Bt);		// could possibly use MAT_REUSE_MATRIX


	MatAXPY(shell->S_approx,-1.0,B_diag_F_inv_Bt,DIFFERENT_NONZERO_PATTERN);		// could possibly use MAT_REUSE_MATRIX
	//MatScale(shell->S_approx,-1.0);		// could possibly use MAT_REUSE_MATRIX


	// ********** SETUP TEMPORARY VECTORS *********************** //
	// shell->temp_vec is u_star
	// shell->temp_vec_2 is delta_p
	// shell->temp_vec_3 is Bt*delta_p
	MatGetVecs(shell->A,&shell->temp_vec,&shell->temp_vec_3);
	MatGetVecs(shell->C,&shell->temp_vec_2,NULL);
	
	

	//ierr = MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);

	// ********* SETUP SCHUR KSP ****************	//
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_schur_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_schur_ksp,"navier_stokes_schur_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_schur_ksp); CHKERRQ(ierr);

	// setup up the operators for the schur solve
	ierr = KSPSetOperators(shell->inner_schur_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);
	
	//ierr = KSPSetType(shell->inner_schur_ksp, KSPPREONLY); CHKERRQ(ierr);
	ierr = KSPSetType(shell->inner_schur_ksp, KSPGMRES); CHKERRQ(ierr);

	PC local_schur_pc;
	ierr = KSPGetPC(shell->inner_schur_ksp,&local_schur_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_schur_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERSUPERLU_DIST);
	//ierr = PCSetType(local_schur_pc,PCHYPRE); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERMUMPS);



	// ************** DESTROY TEMP VECS AND MATS ********************* //
	VecDestroy(&F_diag_inv);
	MatDestroy(&Bt_scaled);
	MatDestroy(&B_diag_F_inv_Bt);
	printf ("done shell pc SIMPLE james setup");

	ierr = PetscLogStagePop();
  return 0;
}





PetscErrorCode SIMPLECShellPCSetUp(PC pc, IS velocity_is, IS pressure_is, KSP schur_ksp)
{
	printf ("inside shell pc SIMPLE james setup");
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);




  ierr = PCShellGetContext(pc,(void**)&shell); CHKERRQ(ierr);


	// ************* WHAT DO WE NEED ****************** //
	// we need a ksp for the velocity solve
	// we need a ksp for the schur solve


	// ************** CREATE THE SUB MATRICES WE WILL USE ********* //
	shell->velocity_is = velocity_is;
	shell->pressure_is = pressure_is;

	Mat A;
	ierr = PCGetOperators(pc,&A,NULL); CHKERRQ(ierr);

  	MatGetSubMatrix(A,shell->velocity_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->A);
  	MatGetSubMatrix(A,shell->velocity_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->Bt);
  	MatGetSubMatrix(A,shell->pressure_is,shell->velocity_is,MAT_INITIAL_MATRIX,&shell->B);
  	MatGetSubMatrix(A,shell->pressure_is,shell->pressure_is,MAT_INITIAL_MATRIX,&shell->C);


	// ********* SETUP VELOCITY KSP ****************** //
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_velocity_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_velocity_ksp,"convection_diffusion_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_velocity_ksp); CHKERRQ(ierr);

	// setup up the operators for the velocity solve
	ierr = KSPSetOperators(shell->inner_velocity_ksp,shell->A,shell->A); CHKERRQ(ierr);
	
	ierr = KSPSetType(shell->inner_velocity_ksp, KSPPREONLY); CHKERRQ(ierr);

	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_velocity_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_velocity_pc,MATSOLVERSUPERLU_DIST);

	// ************ SETUP SCHUR MATRIX ******************** //
	
	// ********** SETUP VELOCITY  DIAGONAL INVERSE ******************* //
	Vec F_diag_inv;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	//MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	//VecReciprocal(F_diag_inv);

	// multiply by transpose, extract diagonal, square root diagonal


	// instead of the diagonal we use the absolute value of the rowsum, this doesn't work so try sum of abs values
	PetscInt          rstart,rend;
	PetscInt          row,ncols,j,nrows;
	PetscReal         val;
	const PetscInt    *cols;
	const PetscScalar *vals;
	MatGetOwnershipRange(shell->A,&rstart,&rend);
	nrows = 0;
	for (row=rstart; row<rend; row++) {
		MatGetRow(shell->A,row,&ncols,&cols,&vals);
		val  = 0.0;
		for (j=0; j<ncols; j++) {
			val += PetscAbsScalar(vals[j]);
		}
		VecSetValue(F_diag_inv,row,val,ADD_VALUES);

		MatRestoreRow(A,row,&ncols,&cols,&vals);
	}
	VecAssemblyBegin(F_diag_inv);
	VecAssemblyEnd(F_diag_inv);

	//std::cout << "F_diag  sum" << std::endl;
	//VecView(F_diag_inv,PETSC_VIEWER_STDOUT_SELF);

	VecReciprocal(F_diag_inv);

	//std::cout << "F_diag_inv sum" << std::endl;
	//VecView(F_diag_inv,PETSC_VIEWER_STDOUT_SELF);


	//std::cout << "A" << std::endl;
	//MatView(shell->A,PETSC_VIEWER_STDOUT_SELF);



	// ********** SETUP S_approx = C -  B* diag(F)^-1* B^T *************************	//
	Mat Bt_scaled;
	Mat B_diag_F_inv_Bt;
	MatDuplicate(shell->Bt,MAT_COPY_VALUES,&Bt_scaled);	// create temporary place for Bt_scaled
	MatDuplicate(shell->C,MAT_COPY_VALUES,&B_diag_F_inv_Bt);	// create sapce for B_diag_F_inv_Bt
	MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);	// set values of S_approx
	MatDiagonalScale(Bt_scaled,F_diag_inv,NULL);

	MatMatMult(shell->B,Bt_scaled,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&B_diag_F_inv_Bt);		// could possibly use MAT_REUSE_MATRIX


	MatAXPY(shell->S_approx,-1.0,B_diag_F_inv_Bt,DIFFERENT_NONZERO_PATTERN);		// could possibly use MAT_REUSE_MATRIX
	//MatScale(shell->S_approx,-1.0);		// could possibly use MAT_REUSE_MATRIX


	// ********** SETUP TEMPORARY VECTORS *********************** //
	// shell->temp_vec is u_star
	// shell->temp_vec_2 is delta_p
	// shell->temp_vec_3 is Bt*delta_p
	MatGetVecs(shell->A,&shell->temp_vec,&shell->temp_vec_3);
	MatGetVecs(shell->C,&shell->temp_vec_2,NULL);
	
	

	//ierr = MatDuplicate(shell->C,MAT_COPY_VALUES,&shell->S_approx);

	// ********* SETUP SCHUR KSP ****************	//
	ierr = KSPCreate(PETSC_COMM_WORLD,&shell->inner_schur_ksp); CHKERRQ(ierr);

	ierr = KSPSetOptionsPrefix(shell->inner_schur_ksp,"navier_stokes_schur_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (shell->inner_schur_ksp); CHKERRQ(ierr);

	// setup up the operators for the schur solve
	ierr = KSPSetOperators(shell->inner_schur_ksp,shell->S_approx,shell->S_approx); CHKERRQ(ierr);
	
	//ierr = KSPSetType(shell->inner_schur_ksp, KSPPREONLY); CHKERRQ(ierr);
	ierr = KSPSetType(shell->inner_schur_ksp, KSPGMRES); CHKERRQ(ierr);

	PC local_schur_pc;
	ierr = KSPGetPC(shell->inner_schur_ksp,&local_schur_pc); CHKERRQ(ierr);

	ierr = PCSetType(local_schur_pc,PCLU); CHKERRQ(ierr);
	ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERSUPERLU_DIST);
	//ierr = PCSetType(local_schur_pc,PCHYPRE); CHKERRQ(ierr);
	//ierr = PCFactorSetMatSolverPackage(local_schur_pc,MATSOLVERMUMPS);



	// ************** DESTROY TEMP VECS AND MATS ********************* //
	VecDestroy(&F_diag_inv);
	MatDestroy(&Bt_scaled);
	MatDestroy(&B_diag_F_inv_Bt);
	printf ("done shell pc SIMPLE james setup");

	ierr = PetscLogStagePop();
  return 0;
}












/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleShellPCApply"
/*
   SampleShellPCApply - This routine demonstrates the use of a
   user-provided preconditioner.

   Input Parameters:
+  pc - preconditioner object
-  x - input vector

   Output Parameter:
.  y - preconditioned vector

   Notes:
   This code implements the Jacobi preconditioner, merely as an
   example of working with a PCSHELL.  Note that the Jacobi method
   is already provided within PETSc.
*/

PetscErrorCode PressureShellPCApply(PC pc,Vec x,Vec y)
{

  NSShellPC  *shell;
  PetscErrorCode ierr;

	ierr = PetscLogStagePush(1);

	//printf ("in Pressure preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	//VecView(x,PETSC_VIEWER_STDOUT_SELF);

	// setup up the operators for the first inner solve

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);
	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);


	ierr = PCApply(local_mass_pc,x,y); CHKERRQ(ierr);
	PetscInt num_mass_its = 0;
  ierr = KSPGetIterationNumber (local_mass_ksp, &num_mass_its);
	//printf ("num_mass_its = %d\n",num_mass_its);


	//ierr = MatMult(shell->pressure_mass_matrix,x,y);

	//VecView(y,PETSC_VIEWER_STDOUT_SELF);

	ierr = PetscLogStagePop();
  return 0;
}


PetscErrorCode PCDShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	PetscLogStage pcd_stage;
	ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	//MatView(pressure_mass_matrix,PETSC_VIEWER_STDOUT_SELF);

	// get the operators for the first inner solve

	PC local_lap_pc;
	ierr = KSPGetPC(shell->inner_lap_ksp,&local_lap_pc); CHKERRQ(ierr);
	KSP local_lap_ksp;
	ierr = PCKSPGetKSP(local_lap_pc,&local_lap_ksp); CHKERRQ(ierr);
	
	ierr = PCApply(local_lap_pc,x,shell->temp_vec);CHKERRQ(ierr);
	PetscInt num_lap_its = 0;
  ierr = KSPGetIterationNumber (local_lap_ksp, &num_lap_its); CHKERRQ(ierr);
	//std::cout << "num_lap_its = " << num_lap_its << std::endl;
	//printf ("num_lap_its = %d\n",num_lap_its);

	PetscInt rows;
	PetscInt cols;
	ierr = MatGetSize(shell->pressure_laplacian_matrix,&rows,&cols);




	// ******* Apply F_p ************ //
	//
	ierr = MatMult(shell->pressure_convection_diffusion_matrix,shell->temp_vec,shell->temp_vec_2); CHKERRQ(ierr);
	//ierr = MatMult(shell->pressure_convection_diffusion_matrix,x,shell->temp_vec_2); CHKERRQ(ierr);

	
	

	// ******* Apply M_p^-1 ********* //

	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);
	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);

	ierr = PCApply(local_mass_pc,shell->temp_vec_2,y); CHKERRQ(ierr);
	PetscInt num_mass_its = 0;
  ierr = KSPGetIterationNumber (local_mass_ksp, &num_mass_its); CHKERRQ(ierr);
	//std::cout << "num_mass_its = " << num_mass_its << std::endl;
	//printf ("num_mass_its = %d\n",num_mass_its);

	ierr = PetscLogStagePop();

  return 0;
}

PetscErrorCode PCD2ShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	PetscLogStage pcd2_stage;
	ierr = PetscLogStagePush(1);

	printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	//MatView(pressure_mass_matrix,PETSC_VIEWER_STDOUT_SELF);

	// setup up the operators for the first inner solve
	PC local_mass_pc;
	ierr = KSPGetPC(shell->inner_mass_ksp,&local_mass_pc); CHKERRQ(ierr);

	// set the pc to be a ksp
	// should probably use gmres by default okay
	ierr = PCSetType(local_mass_pc,PCKSP); CHKERRQ(ierr);


	KSP local_mass_ksp;
	ierr = PCKSPGetKSP(local_mass_pc,&local_mass_ksp); CHKERRQ(ierr);


	ierr = KSPSetOptionsPrefix(local_mass_ksp,"mass_matrix_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_mass_ksp); CHKERRQ(ierr);
	

	ierr = KSPSetOperators(local_mass_ksp,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);
	ierr = PCSetOperators(local_mass_pc,shell->pressure_mass_matrix,shell->pressure_mass_matrix); CHKERRQ(ierr);

	
	ierr = PCApply(local_mass_pc,x,shell->temp_vec);CHKERRQ(ierr);
	PetscInt num_mass_its = 0;
  ierr = KSPGetIterationNumber (local_mass_ksp, &num_mass_its); CHKERRQ(ierr);
	std::cout << "num_mass_its = " << num_mass_its << std::endl;

	//ierr = MatMult(shell->pressure_laplacian_matrix,x,shell->temp_vec_2); CHKERRQ(ierr);
	







	// ******* Apply F_p ************ //
	//
	ierr = MatMult(shell->pressure_convection_diffusion_matrix,shell->temp_vec,shell->temp_vec_2); CHKERRQ(ierr);
	//ierr = MatMult(shell->pressure_convection_diffusion_matrix,x,shell->temp_vec_2); CHKERRQ(ierr);

	
	

	// ******* Apply A_p^-1 ********* //

	PC local_lap_pc;
	ierr = KSPGetPC(shell->inner_lap_ksp,&local_lap_pc); CHKERRQ(ierr);
	ierr = PCSetType(local_lap_pc,PCKSP); CHKERRQ(ierr);
	KSP local_lap_ksp;
	ierr = PCKSPGetKSP(local_lap_pc,&local_lap_ksp); CHKERRQ(ierr);

	
	ierr = KSPSetOptionsPrefix(local_lap_ksp,"laplacian_"); CHKERRQ(ierr);
	ierr = KSPSetFromOptions (local_lap_ksp); CHKERRQ(ierr);
	

	ierr = KSPSetOperators(local_lap_ksp,shell->pressure_laplacian_matrix,shell->pressure_laplacian_preconditioner); CHKERRQ(ierr);
	ierr = PCSetOperators(local_lap_pc,shell->pressure_laplacian_matrix,shell->pressure_laplacian_preconditioner); CHKERRQ(ierr);

	ierr = PCApply(local_lap_pc,shell->temp_vec_2,y); CHKERRQ(ierr);
	//PCApply(local_lap_pc,x,y);
	PetscInt num_lap_its = 0;
  ierr = KSPGetIterationNumber (local_lap_ksp, &num_lap_its); CHKERRQ(ierr);
	std::cout << "num_lap_its = " << num_lap_its << std::endl;
	PetscReal rnorm = 0.;
	KSPGetResidualNorm(local_lap_ksp,&rnorm);
	std::cout << "norm = " << rnorm << std::endl;
	KSPConvergedReason creason;
	KSPGetConvergedReason(local_lap_ksp,&creason);
	std::cout << "creason = " << creason << std::endl;
	
	ierr = PetscLogStagePop();

  return 0;
}


PetscErrorCode MonolithicShellPCApply(PC pc,Vec x,Vec y)
{

  	NSShellPC  *shell;
  	PetscErrorCode ierr;

	//ierr = PetscLogStagePush(2);


  	ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	
	// apply the preconditioner (lu)	
	PC local_velocity_pc;
	ierr = KSPGetPC(shell->inner_velocity_ksp,&local_velocity_pc); CHKERRQ(ierr);
	ierr = PCApply(local_velocity_pc,x,y); CHKERRQ(ierr);




	//ierr = PetscLogStagePop();
  return 0;
}










PetscErrorCode LSCShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);

	Mat A,B,C;
	MatSchurComplementGetSubMatrices(pcmat,&A,NULL,&B,&C,NULL);

	KSPSolve(shell->inner_lap_ksp,x,shell->temp_vec_3);

	MatMult(B,shell->temp_vec_3,shell->temp_vec);

	//VecPointwiseMult(shell->temp_vec,shell->temp_vec,shell->lsc_scale);

	MatMult(A,shell->temp_vec,shell->temp_vec_2);

	//VecPointwiseMult(shell->temp_vec_2,shell->temp_vec_2,shell->lsc_scale);

	MatMult(C,shell->temp_vec_2,shell->temp_vec_3);
	KSPSolve(shell->inner_lap_ksp,shell->temp_vec_3,y);
	return(0);



	ierr = PetscLogStagePop();

  return 0;
}






PetscErrorCode LSCScaledShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);

	Mat A,B,C;
	MatSchurComplementGetSubMatrices(pcmat,&A,NULL,&B,&C,NULL);

	KSPSolve(shell->inner_lap_ksp,x,shell->temp_vec_3);

	MatMult(B,shell->temp_vec_3,shell->temp_vec);

	VecPointwiseMult(shell->temp_vec,shell->temp_vec,shell->lsc_scale);

	MatMult(A,shell->temp_vec,shell->temp_vec_2);

	VecPointwiseMult(shell->temp_vec_2,shell->temp_vec_2,shell->lsc_scale);

	MatMult(C,shell->temp_vec_2,shell->temp_vec_3);
	KSPSolve(shell->inner_lap_ksp,shell->temp_vec_3,y);
	return(0);



	ierr = PetscLogStagePop();

  return 0;
}





PetscErrorCode LSCScaledStabilisedShellPCApply(PC pc,Vec x,Vec y)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	//printf ("inside shell pc lsc scaled stabilised james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);

	Mat A,B,C;
	MatSchurComplementGetSubMatrices(pcmat,&A,NULL,&B,&C,NULL);

	KSPSolve(shell->inner_lap_ksp,x,shell->temp_vec_3);

	//std::cout << "x: " << std::endl;
	double norm;
	VecNorm(x,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	//std::cout << "shell->temp_vec_3: " << std::endl;
	VecNorm(shell->temp_vec_3,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatMult(B,shell->temp_vec_3,shell->temp_vec);

	//std::cout << "shell->temp_vec: " << std::endl;
	VecNorm(shell->temp_vec,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecPointwiseMult(shell->temp_vec,shell->temp_vec,shell->lsc_scale);

	//std::cout << "shell->temp_vec: " << std::endl;
	VecNorm(shell->temp_vec,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	MatMult(A,shell->temp_vec,shell->temp_vec_2);

	//std::cout << "shell->temp_vec_2: " << std::endl;
	VecNorm(shell->temp_vec_2,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecPointwiseMult(shell->temp_vec_2,shell->temp_vec_2,shell->lsc_scale);

	//std::cout << "shell->temp_vec_2: " << std::endl;
	VecNorm(shell->temp_vec_2,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	MatMult(C,shell->temp_vec_2,shell->temp_vec_3);

	//std::cout << "shell->temp_vec_3: " << std::endl;
	VecNorm(shell->temp_vec_3,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	KSPSolve(shell->inner_lap_ksp,shell->temp_vec_3,y);

	//std::cout << "y: " << std::endl;
	VecNorm(y,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

//lsc_stab_alpha_D_inv
	// + alpha D^-1
	//Vec alpha_D_inv_x;
	VecPointwiseMult(shell->temp_vec_3,shell->lsc_stab_alpha_D_inv,x);

	//std::cout << "shell->temp_vec_3: " << std::endl;
	VecNorm(shell->temp_vec_3,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecAXPY(y,1.0,shell->temp_vec_3);
	
	//std::cout << "y: " << std::endl;
	VecNorm(y,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}




PetscErrorCode SIMPLEShellPCApply(PC pc,Vec x,Vec y)
{
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	// u_star is u_star
	// delta_p is delta_p
	// F_diag_inv_Bt_delta_p is Bt*delta_p
	// r_u is r_u
	// r_p is r_p
	// B_u_star is B*u_star
	// y_u is y_u
	// y_p is y_p

	// get some matrices we will need
	// setup some temporary vectors
	Vec F_diag_inv,delta_p,F_diag_inv_Bt_delta_p,r_u,r_p,y_u,y_p,r_p_minus_B_u_star,u_star;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	MatGetVecs(shell->A,&u_star,NULL);
	MatGetVecs(shell->C,&delta_p,NULL);
	MatGetVecs(shell->A,&F_diag_inv_Bt_delta_p,NULL);
	MatGetVecs(shell->C,&r_p_minus_B_u_star,NULL);


	// don't need to allocate
	//MatGetVecs(shell->A,&r_u,NULL);
	//MatGetVecs(shell->C,&r_p,NULL);
	VecGetSubVector(x,shell->velocity_is,&r_u);
	VecGetSubVector(x,shell->pressure_is,&r_p);

	//MatGetVecs(shell->A,&y_u,NULL);
	//MatGetVecs(shell->C,&y_p,NULL);
	VecGetSubVector(y,shell->velocity_is,&y_u);
	VecGetSubVector(y,shell->pressure_is,&y_p);

	//std::cout << "r_u:" << std::endl; 
	double norm;
	VecNorm(r_u,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	//std::cout << "r_p:" << std::endl; 
	VecNorm(r_p,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(F_diag_inv);

	//std::cout << "F_diag_inv:" << std::endl; 
	VecNorm(F_diag_inv,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

		

	

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);


	// F^-1 * u_star = r_u
	KSPSolve(shell->inner_velocity_ksp,r_u,u_star);

	//std::cout << "u_star:" << std::endl;
	VecNorm(u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// S_approx * delta_p = r_p - B * u_star
	MatMult(shell->B,u_star,r_p_minus_B_u_star);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecAXPY(r_p_minus_B_u_star,-1.0,r_p);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecScale(r_p_minus_B_u_star,-1.0);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	KSPSolve(shell->inner_schur_ksp,r_p_minus_B_u_star,delta_p);

	//std::cout << "delta_p:" << std::endl; 
	VecNorm(delta_p,NORM_2,&norm); 
	//VecView(delta_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// u = u_star - diag(F)^-1 * Bt * delta_p
	VecCopy(u_star,y_u);
	MatMult(shell->Bt,delta_p,F_diag_inv_Bt_delta_p);
	VecPointwiseMult(F_diag_inv_Bt_delta_p,F_diag_inv,F_diag_inv_Bt_delta_p);
	VecAXPY(y_u,-1.0,F_diag_inv_Bt_delta_p);

	//std::cout << "y_u:" << std::endl; 
	VecNorm(y_u,NORM_2,&norm); 
	//VecView(y_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// p = delta_p
	VecCopy(delta_p,y_p);

	//std::cout << "y_p:" << std::endl; 
	VecNorm(y_p,NORM_2,&norm); 
	//VecView(y_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// put the subvector back in its place
	VecRestoreSubVector(y,shell->velocity_is,&y_u);
	VecRestoreSubVector(y,shell->pressure_is,&y_p);


	//std::cout << "y:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	VecCopy(x,y);

	//std::cout << "y after:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	

	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}





PetscErrorCode SIMPLERShellPCApply(PC pc,Vec x,Vec y)
{
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;

	printf ("inside shell pc lsc james pcapply");

	PetscLogStage pcd_stage;
	//ierr = PetscLogStagePush(1);

	//printf ("in PCD preconditioner\n");

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);

	// u_star is u_star
	// delta_p is delta_p
	// F_diag_inv_Bt_delta_p is Bt*delta_p
	// r_u is r_u
	// r_p is r_p
	// B_u_star is B*u_star
	// y_u is y_u
	// y_p is y_p

	// get some matrices we will need
	// setup some temporary vectors
	Vec F_diag_inv,delta_p,F_diag_inv_Bt_delta_p,r_u,r_p,y_u,y_p,r_p_minus_B_u_star,u_star,r_u_minus_Bt_p_star,p_star,r_p_minus_B_F_diag_inv_r_u,F_diag_inv_r_u,C_p_star;
	MatGetVecs(shell->A,&F_diag_inv,NULL);
	MatGetVecs(shell->A,&u_star,NULL);
	MatGetVecs(shell->C,&delta_p,NULL);
	MatGetVecs(shell->A,&F_diag_inv_Bt_delta_p,NULL);
	MatGetVecs(shell->C,&r_p_minus_B_u_star,NULL);
	MatGetVecs(shell->A,&r_u_minus_Bt_p_star,NULL);
	MatGetVecs(shell->C,&r_p_minus_B_F_diag_inv_r_u,NULL);
	MatGetVecs(shell->A,&F_diag_inv_r_u,NULL);
	MatGetVecs(shell->C,&p_star,NULL);
	MatGetVecs(shell->C,&C_p_star,NULL);


	// don't need to allocate
	//MatGetVecs(shell->A,&r_u,NULL);
	//MatGetVecs(shell->C,&r_p,NULL);
	VecGetSubVector(x,shell->velocity_is,&r_u);
	VecGetSubVector(x,shell->pressure_is,&r_p);

	//MatGetVecs(shell->A,&y_u,NULL);
	//MatGetVecs(shell->C,&y_p,NULL);
	VecGetSubVector(y,shell->velocity_is,&y_u);
	VecGetSubVector(y,shell->pressure_is,&y_p);

	//std::cout << "r_u:" << std::endl; 
	double norm;
	VecNorm(r_u,NORM_2,&norm);
	//VecView(r_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	//std::cout << "r_p:" << std::endl; 
	VecNorm(r_p,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatGetDiagonal(shell->A,F_diag_inv); /* Should be the mass matrix, but we don't have plumbing for that yet */
	VecReciprocal(F_diag_inv);

	//std::cout << "F_diag_inv:" << std::endl; 
	VecNorm(F_diag_inv,NORM_2,&norm);
	//VecView(r_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

		

	

	Mat pcmat;
	ierr = PCGetOperators(pc,&pcmat,NULL); CHKERRQ(ierr);


	// S_approx * p_star = r_p - B * F_diag^-1* r_u
	VecCopy(r_u,F_diag_inv_r_u);
	VecPointwiseMult(F_diag_inv_r_u,F_diag_inv,F_diag_inv_r_u);
	MatMult(shell->B,F_diag_inv_r_u,r_p_minus_B_F_diag_inv_r_u);
	VecAXPY(r_p_minus_B_F_diag_inv_r_u,-1.0,r_p);
	VecScale(r_p_minus_B_F_diag_inv_r_u,-1.0);
	KSPSolve(shell->inner_schur_ksp,r_p_minus_B_F_diag_inv_r_u,p_star);
	
	//std::cout << "p_star:" << std::endl;
	VecNorm(p_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	//std::cout << "r_p_minus_B_F_diag_inv_r_u:" << std::endl;
	VecNorm(r_p_minus_B_F_diag_inv_r_u,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	

	// F^-1 * u_star = r_u - Bt * p_star
	MatMult(shell->Bt,p_star,r_u_minus_Bt_p_star);
	VecAXPY(r_u_minus_Bt_p_star,-1.0,r_u);
	VecScale(r_u_minus_Bt_p_star,-1.0);
	
	KSPSolve(shell->inner_velocity_ksp,r_u_minus_Bt_p_star,u_star);

	//std::cout << "u_star:" << std::endl;
	VecNorm(u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// S_approx * delta_p = r_p - B * u_star - C*p_star
	MatMult(shell->B,u_star,r_p_minus_B_u_star);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecAXPY(r_p_minus_B_u_star,-1.0,r_p);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;
	VecScale(r_p_minus_B_u_star,-1.0);
	//std::cout << "r_p_minus_B_u_star:" << std::endl;
	VecNorm(r_p_minus_B_u_star,NORM_2,&norm); 
	//VecView(u_star,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	MatMult(shell->C,p_star,C_p_star);
	VecAXPY(r_p_minus_B_u_star,-1.0,C_p_star);
	KSPSolve(shell->inner_schur_ksp,r_p_minus_B_u_star,delta_p);

	//std::cout << "delta_p:" << std::endl; 
	VecNorm(delta_p,NORM_2,&norm); 
	//VecView(delta_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// u = u_star - diag(F)^-1 * Bt * delta_p
	VecCopy(u_star,y_u);
	MatMult(shell->Bt,delta_p,F_diag_inv_Bt_delta_p);
	VecPointwiseMult(F_diag_inv_Bt_delta_p,F_diag_inv,F_diag_inv_Bt_delta_p);
	VecAXPY(y_u,-1.0,F_diag_inv_Bt_delta_p);

	//std::cout << "y_u:" << std::endl; 
	VecNorm(y_u,NORM_2,&norm); 
	//VecView(y_u,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// p = delta_p
	VecCopy(delta_p,y_p);
	VecAXPY(y_p,1.0,p_star);
	

	//std::cout << "y_p:" << std::endl; 
	VecNorm(y_p,NORM_2,&norm); 
	//VecView(y_p,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	// put the subvector back in its place
	VecRestoreSubVector(y,shell->velocity_is,&y_u);
	VecRestoreSubVector(y,shell->pressure_is,&y_p);


	//std::cout << "y:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	//VecCopy(x,y);

	//std::cout << "y after:" << std::endl; 
	VecNorm(y,NORM_2,&norm); 
	//VecView(y,PETSC_VIEWER_STDOUT_SELF);
	//std::cout << norm << std::endl;

	
	

	return(0);



	//ierr = PetscLogStagePop();

  return 0;
}


/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleNSShellDestroy"
/*
   SampleNSShellDestroy - This routine destroys a user-defined
   preconditioner context.

   Input Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode NSShellDestroy(PC pc)
{
  NSShellPC  *shell;
  PetscErrorCode ierr;

	//ierr = PetscLogStagePush(1);
  std::cout << "in shell destroy" << std::endl;

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);
	ierr = KSPDestroy(&shell->inner_mass_ksp);CHKERRQ(ierr);
	ierr = KSPDestroy(&shell->inner_lap_ksp);CHKERRQ(ierr);
	ierr = KSPDestroy(&shell->inner_velocity_ksp);CHKERRQ(ierr);
	ierr = KSPDestroy(&shell->inner_schur_ksp);CHKERRQ(ierr);
	ierr = VecDestroy(&shell->temp_vec);CHKERRQ(ierr);
	ierr = VecDestroy(&shell->temp_vec_2);CHKERRQ(ierr);
	ierr = VecDestroy(&shell->temp_vec_3);CHKERRQ(ierr);
	ierr = VecDestroy(&shell->lsc_scale);CHKERRQ(ierr);
	ierr = MatDestroy(&shell->S_approx);CHKERRQ(ierr);
//	ierr = MatDestroy(&shell->velocity_matrix);CHKERRQ(ierr);
	ierr = MatDestroy(&shell->lsc_laplacian_matrix);CHKERRQ(ierr);
	ierr = VecDestroy(&shell->lsc_stab_alpha_D_inv);CHKERRQ(ierr);


  ierr = PetscFree(shell);CHKERRQ(ierr);
	
  std::cout << "done shell destroy" << std::endl;
	//ierr = PetscLogStagePop();

  return 0;
}




/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "SampleNSShellDestroy"
/*
   SIMPLEShellDestroy - This routine destroys a user-defined
   preconditioner context.

   Input Parameter:
.  shell - user-defined preconditioner context
*/
PetscErrorCode SIMPLEShellDestroy(PC pc)
{
  SIMPLEShellPC  *shell;
  PetscErrorCode ierr;


	//ierr = PetscLogStagePush(1);
  std::cout << "in SIMPLE shell destroy" << std::endl;

  ierr = PCShellGetContext(pc,(void**)&shell);CHKERRQ(ierr);
	ierr = KSPDestroy(&shell->inner_schur_ksp);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = KSPDestroy(&shell->inner_velocity_ksp);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec_2);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = VecDestroy(&shell->temp_vec_3);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->S_approx);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl; //here, it is getting deleted earlier...
//	ierr = MatDestroy(&shell->velocity_matrix);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->A);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->Bt);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->B);CHKERRQ(ierr);
  std::cout << "in SIMPLE shell destroy" << std::endl;
	ierr = MatDestroy(&shell->C);CHKERRQ(ierr);



  ierr = PetscFree(shell);CHKERRQ(ierr);
	
  std::cout << "done SIMPLE shell destroy" << std::endl;
	//ierr = PetscLogStagePop();

  return 0;
}



// ************ MATSHELL MULT ***************** //
PetscErrorCode MatShellMultFull(Mat A, Vec vx, Vec vy)
{
	PetscErrorCode ierr;
	PCD2ShellMatrixCtx *ctx;


	ierr = MatShellGetContext(A, (void **)&ctx); CHKERRQ(ierr);

	KSP schur_ksp = ctx->schur_ksp;

	// setup up the operators for the first inner solve
	Mat S, Pmat;
	//MatStructure mat_structure;

	//MatView(pressure_laplacian_matrix->mat(),PETSC_VIEWER_STDOUT_SELF);
	//ierr = KSPGetOperators(schur_ksp,&S,NULL,NULL); CHKERRQ(ierr);
	ierr = KSPGetOperators(schur_ksp,&S,NULL); CHKERRQ(ierr);

	Mat A00;
	Mat Bt;
	Mat B;

	//ierr = MatSchurComplementGetSubmatrices(S,&A00,NULL,&Bt,&B,NULL);
	ierr = MatSchurComplementGetSubMatrices(S,&A00,NULL,&Bt,&B,NULL);

	//MatView(S,PETSC_VIEWER_STDOUT_SELF);
	//MatView(B,PETSC_VIEWER_STDOUT_SELF);

	//std::cout << "BT baby" << std::endl;
	//MatView(Bt,PETSC_VIEWER_STDOUT_SELF);

	//std::cout << "not made vecs" << std::endl;
	Vec temp_vec_1;
	Vec temp_vec_2;
	Vec inv_diag;
  // create some vectors for temporary pc applying
	PetscInt n_row;
	PetscInt n_row_local;
	//std::cout << "not got mat" << std::endl;
	ierr = MatGetSize(Bt,&n_row,NULL); CHKERRQ(ierr);
	ierr = MatGetLocalSize(Bt,&n_row_local,NULL); CHKERRQ(ierr);
	//std::cout << "got mat sizes rows: " << n_row << ", cols = " << n_col << std::endl;
	// create some vectors that can hold temporary info
	ierr = VecCreate(PETSC_COMM_WORLD,&temp_vec_1); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&temp_vec_2); CHKERRQ(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&inv_diag); CHKERRQ(ierr);
	ierr = VecSetSizes(temp_vec_1, n_row_local, n_row); CHKERRQ(ierr);
	ierr = VecSetSizes(temp_vec_2, n_row_local, n_row); CHKERRQ(ierr);
	ierr = VecSetSizes(inv_diag, n_row_local, n_row); CHKERRQ(ierr);
	ierr = VecSetFromOptions(temp_vec_1); CHKERRQ(ierr);
	ierr = VecSetFromOptions(temp_vec_2); CHKERRQ(ierr);
	ierr = VecSetFromOptions(inv_diag); CHKERRQ(ierr);

	ierr = MatGetDiagonal(ctx->velocity_mass_matrix,inv_diag);CHKERRQ(ierr);
	//VecView(inv_diag,PETSC_VIEWER_STDOUT_SELF);
	ierr = VecReciprocal(inv_diag);CHKERRQ(ierr);

	PetscInt size;
	ierr = VecGetSize(temp_vec_1, &size); CHKERRQ(ierr);
	//std::cout << "size1 =  " << size << std::endl; 
	ierr = VecGetSize(temp_vec_2, &size); CHKERRQ(ierr);
	//std::cout << "size2 =  " << size << std::endl; 

	//ierr = MatGetSize(ctx->velocity_mass_matrix,&n_row,&n_col); CHKERRQ(ierr);
	//std::cout << "got mat sizes rows: " << n_row << ", cols = " << n_col << std::endl;

	//VecView(inv_diag,PETSC_VIEWER_STDOUT_SELF);
	//MatView(Bt,PETSC_VIEWER_STDOUT_SELF);
	//MatView(ctx->velocity_mass_matrix,PETSC_VIEWER_STDOUT_SELF);
	//MatView(B,PETSC_VIEWER_STDOUT_SELF);

	//std::cout << "made vecs" << std::endl;
	ierr = MatMult(Bt,vx,temp_vec_1); CHKERRQ(ierr);
	//std::cout << "1" << std::endl;
	//ierr = MatMult(A00,temp_vec_1,temp_vec_2); CHKERRQ(ierr);
	ierr = VecPointwiseMult(temp_vec_2,inv_diag,temp_vec_1);
	//std::cout << "2" << std::endl;
	ierr = MatMult(B,temp_vec_2,vy); CHKERRQ(ierr);
	//std::cout << "3" << std::endl;

	ierr = VecDestroy(&temp_vec_1);CHKERRQ(ierr);
	ierr = VecDestroy(&temp_vec_2);CHKERRQ(ierr);
	ierr = VecDestroy(&inv_diag);CHKERRQ(ierr);

	//ierr = MatMult(ctx->pressure_laplacian_matrix,vx,vy); CHKERRQ(ierr);

	return 0;
}

