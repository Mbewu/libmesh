#!/bin/bash

# this is a template file that is used to generate a test_script that either runs a file or submits a job
# #######
# parameters to be passed in
# - COMPUTER_TYPE - the type of computer you are running on (laptop, compute-lung, arcus-a or arcus-b)
# - NUM_PROCS - number of procs, for non-batch computer
# - QUEUE_TYPE - type of queue, devel for test queue, empty string for normal queue
# - JOB_NAME - name of the job
# - WALLTIME - max wall time in format "100:00:00"
# - NUM_NODES - number of nodes on a batch system
# - NUM_PROCS_PER_NODE - number of procs per node on a batch system
# - OUTPUT_DIR - the directory for the results to go in
# - LIBMESH_OPTIONS - the options to put to libmesh

##################### USER PARAMETERS ###########################
# ideally these parameters should be read in as parameters and then the script can do its work

COMPUTER_TYPE="laptop"      #"laptop" is laptop, "compute-lung" is compute-lung, "arcus-a" is arcus-a and "arcus-b" is arcus-b

# non-batch specific variables
NUM_PROCS=1

# batch specific variables
QUEUE_TYPE="devel"           #"devel" is test and "" is just normal
JOB_NAME="arc_job"
WALLTIME=00:10:00
NUM_NODES=1
NUM_PROCS_PER_NODE=1

######################################################################################################






########################## SETUP SOME VARIABLES AND PARAMETERS (AUTOMATIC)  ########################

# set the code BASE_DIR and the OUTPUT_BASE_DIR
if [ "$COMPUTER_TYPE" = "laptop" ] ; then
	BASE_DIR="/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes"
elif [ "$COMPUTER_TYPE" = "compute-lung" ] ; then
	BASE_DIR="/users/jmbewu/coupled_navier_stokes"
elif [ "$COMPUTER_TYPE" = "arcus-a" ] ; then
	BASE_DIR="/home/comp-respiratory-modelling/jmbewu/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes"
elif [ "$COMPUTER_TYPE" = "arcus-b" ] ; then
	BASE_DIR="/home/comp-respiratory-modelling/jmbewu/libmesh-git-b/libmesh/examples/dphil/coupled_navier_stokes"
else
	echo "ERROR: invalid COMPUTER_TYPE specified."
fi

# variables for script
THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
THIS_FILE="$0"



#############################################################################





############################ SIMULATION PARAMETERS ###########################

OUTPUT_DIR="/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_20/particle_size_2/"
LIBMESH_OPTIONS="-input_file /home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/libmesh_scripts/particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_20/navier.in -particle_deposition_input_file /home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/libmesh_scripts/particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_20/particle_deposition.in -particle_diameter 2.e-6 -output_folder /home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_20/particle_size_2// --use-petsc --solver_variable_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_w 0 --solver_group_ns3d_p 1  --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_p 0 --solver_group_ns3d1d_w 0 --solver_group_ns3d1d_Q 1 --solver_group_ns3d1d_P 1 --solver_group_ns3d1d_0_u 0 --solver_group_ns3d1d_0_v 0 --solver_group_ns3d1d_0_w 0 --solver_group_ns3d1d_0_p 1 --solver_system_names -log_summary -preconditioner_type_3d 0 -preconditioner_type_3d1d 0 -geometry_type 4 -mesh_file /home/james/meshes-git/dphil_meshes/meshes_for_deposition/multi_bifurcation_symmetric_diff_wall_bdy_2.msh -mesh_input_scaling_3d 0.02 -inflow_bdy_id 1011 -gmsh_diff_wall_bdy_id true -_read_1d_mesh false -num_generations_string 0 -num_generations 8 -twod_oned_tree true -threed false -stokes false -unsteady 1 -sim_type 3 -match_1d_to_3d_mesh true -reynolds_number_calculation false -prescribed_flow 0 -period 4.0 -end_time 2.0 -dt 0.01 -write_interval 0.01 -nonlinear_tolerance 1e-6 -density 1.2 -viscosity 1.8e-5 -newton 1 -bifurcation_start_1d true -half_initial_length false -use_centreline_data true -input_1d_file  /home/james/meshes-git/dphil_meshes/meshes_for_deposition/multi_bifurcation_symmetric_diff_wall_bdy_2 -radius_on_edge false -particle_deposition 6 -unsteady_from_steady false -restart_folder /home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_20/fluid_sim/ -velocity_mag_3d 2.122 -particle_end_time 2.0 -particle_deposition_rate 10000 -particle_dt 0.00005"

###########################################################################








####################### RUN PROGRAM ON NORMAL PC ##########################

if [ "$COMPUTER_TYPE" = "laptop" ] || [ "$COMPUTER_TYPE" = "compute-lung" ] ; then
	$BASE_DIR/run_program.sh "$BASE_DIR" "$LIBMESH_OPTIONS" "$OUTPUT_DIR" "$NUM_PROCS"
	#$BASE_DIR/example-opt $LIBMESH_OPTIONS 2>&1 | tee $OUTPUT_DIR/output.log
fi

###########################################################################




######################## RUN PROGRAM ON arcus-a ##########################

# - create the batch file (need to pass it some variables)
# - submit the batch file and return the job name/number

if [ "$COMPUTER_TYPE" = "arcus-a" ]; then
	# generate file and make executable
	$BASE_DIR/generate_pbs_file.sh "$JOB_NAME" "$WALLTIME" "$NUM_NODES" "$NUM_PROCS_PER_NODE" "$LIBMESH_OPTIONS" "$THIS_DIR" "$BASE_DIR"
	chmod +x $THIS_DIR/job_script.sh

	# submit job and record the job name/number
	if [ "$QUEUE_TYPE" = "devel" ]; then
		JOB_ID=$(qsub -q develq $THIS_DIR/job_script.sh)
	else
		JOB_ID=$(qsub $THIS_DIR/job_script.sh)
	fi
	
 

	# copy the batch file to the output directory
	cp $THIS_DIR/job_script.sh $OUTPUT_DIR/
	# copy the job id to the output directory
	echo "$JOB_ID" > $OUTPUT_DIR/job_id.dat
fi

#####################################################################



######################## RUN PROGRAM ON arcus-b ##########################

# - create the batch file (need to pass it some variables)
# - submit the batch file and return the job name/number

if [ "$COMPUTER_TYPE" = "arcus-b" ]; then
	# generate file and make executable
	$BASE_DIR/generate_slurm_file.sh "$JOB_NAME" "$WALLTIME" "$NUM_NODES" "$NUM_PROCS_PER_NODE" "$LIBMESH_OPTIONS" "$THIS_DIR" "$BASE_DIR"
	chmod +x $THIS_DIR/job_script.sh

	# submit job and record the job name/number (doesn't record job name/number on arcus-b)
	if [ "$QUEUE_TYPE" = "devel" ]; then
		JOB_ID=$(sbatch -p devel $THIS_DIR/job_script.sh)
	else
		JOB_ID=$(sbatch $THIS_DIR/job_script.sh)
	fi
	
 

	# copy the batch file to the output directory
	cp $THIS_DIR/job_script.sh $OUTPUT_DIR/
	# copy the job id to the output directory
	echo "$JOB_ID" > $OUTPUT_DIR/job_id.dat
fi

#####################################################################








################# OUTPUT SOME HOUSEKEEPING DATA #####################

# copy this script to the output folder
cp $THIS_FILE $OUTPUT_DIR/

# output the git version in a file
GIT_VERSION=$(git describe)
echo "$GIT_VERSION" > $OUTPUT_DIR/git_version.dat

# output the name of the computer and the directory it was run from in a file
COMPUTER_NAME=$(hostname)
echo "$COMPUTER_NAME" > $OUTPUT_DIR/computer_name.dat
echo "$THIS_DIR" > $OUTPUT_DIR/run_directory.dat

##############################################################

