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

# Pass environment to compute nodes
#PBS -V

# Give the job a name
#PBS -N %JOB_NAME%

# Walltime limit
#PBS -l walltime=%WALLTIME%

# Use 2 nodes, each running 8 processes (Hal and Sal) or 16 (Arcus)
#PBS -l nodes=%NUM_NODES%:ppn=%NUM_PROCS_PER_NODE%

# Email me at the beginning and end of job
#PBS -m be

# at the following address
#PBS -M james.mbewu@gmail.com

# Go to the directory where you submitted the job
cd $PBS_O_WORKDIR

### script to define MPI_HOSTS for mpirun -- uncomment as appropriate
## Hal and Sal: source the script enable_hal_mpi.sh
##. enable_hal_mpi.sh
## Arcus: source the script enable_arcus_mpi.sh
. enable_arcus_mpi.sh

# Launch the application
/system/software/arcus/mpi/openmpi/1.8.1/gcc-4.8.2/bin/mpirun $MPI_HOSTS ./example-opt %LIBMESH_OPTIONS%
# mpirun -np $MPI_NPROCS -machinefile $PBS_NODEFILE ./cluster_myprog

