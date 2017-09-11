#!/bin/sh

##############################################
## Template for creating PBS files on arcus-a
## - parameters:
## - JOB_NAME - name of job (text)
## - WALLTIME - max walltime (00:00:00)
## - NUM_NODES - number of nodes
## - NUM_PROCS_PER_NODE - number of procs per node
## - LIBMESH_OPTIONS - command line options for libmesh program
##############################################

# Give the job a name
#SBATCH --job-name=arc_job

# Walltime limit
#SBATCH --time=00:10:00

# Use 2 nodes, each running 8 processes (Hal and Sal) or 16 (Arcus)
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16


### script to define MPI_HOSTS for mpirun -- uncomment as appropriate
## Hal and Sal: source the script enable_hal_mpi.sh
##. enable_hal_mpi.sh
## Arcus: source the script enable_arcus_mpi.sh
. enable_arcus-b_mpi.sh

# Launch the application
/system/software/arcus/mpi/openmpi/1.8.1/gcc-4.8.2/bin/mpirun $MPI_HOSTS ./example-opt -input_file /home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/libmesh_scripts/test_script/navier.in -output_folder /home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes//results/reproducibility_testing/test_1/ -sim_type 5 -dt 5.0  --use-petsc --solver_variable_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_w 0 --solver_group_ns3d_p 1  --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_p 0 --solver_group_ns3d1d_w 0 --solver_group_ns3d1d_Q 1 --solver_group_ns3d1d_P 1 --solver_group_ns3d1d_0_u 0 --solver_group_ns3d1d_0_v 0 --solver_group_ns3d1d_0_w 0 --solver_group_ns3d1d_0_p 1 --solver_system_names -log_summary -newton 0 -streamline_diffusion true -stokes true -fieldsplit true -no_refinement 0 -output_no_refinement 0 -linear_shape_functions true -stab true -preconditioner_type_3d1d 12 -preconditioner_type_3d 9 -reuse_preconditioner false -reynolds_number 500 -period 1000.0 -end_time 500.0 -mesh_file /home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes//meshes/meshes_for_paper/multi_bifurcation_symmetric_unstructured_2.msh -num_generations_string 99 -unsteady 0
# mpirun -np $MPI_NPROCS -machinefile $PBS_NODEFILE ./cluster_myprog

