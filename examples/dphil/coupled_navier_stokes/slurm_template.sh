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
#SBATCH --job-name=%JOB_NAME%

# Walltime limit
#SBATCH --time=%WALLTIME%

# Use 2 nodes, each running 8 processes (Hal and Sal) or 16 (Arcus)
#SBATCH --nodes=%NUM_NODES%
#SBATCH --ntasks-per-node=%NUM_PROCS_PER_NODE%


### script to define MPI_HOSTS for mpirun -- uncomment as appropriate
## Hal and Sal: source the script enable_hal_mpi.sh
##. enable_hal_mpi.sh
## Arcus: source the script enable_arcus_mpi.sh
. enable_arcus-b_mpi.sh

# Launch the application
/system/software/arcus-b/lib/mpi/openmpi/1.8.4/gcc-4.9.2/bin/mpirun $MPI_HOSTS %BASE_DIR%/example-opt %LIBMESH_OPTIONS%
# mpirun -np $MPI_NPROCS -machinefile $PBS_NODEFILE ./cluster_myprog

