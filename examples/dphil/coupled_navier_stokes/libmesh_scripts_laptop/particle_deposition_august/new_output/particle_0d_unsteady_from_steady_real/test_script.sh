#!/bin/bash

# in this file we will run a representative simulation
# user must edit USER PARAMETERS and USER SIMULATION PARAMETERS sections

##################### USER PARAMETERS ###########################
# ideally these parameters should be read in as parameters and then the script can do its work

COMPUTER_TYPE="laptop"         #"laptop" is laptop, "compute-lung" is compute-lung, "arcus-a" is arcus-a and "arcus-b" is arcus-b

# non-batch specific variables
NUM_PROCS="1"

# batch specific variables
QUEUE_TYPE="devel"              #"devel" is test and "" is just normal
JOB_NAME="arc_job"
WALLTIME="00:10:00"
NUM_NODES=1
NUM_PROCS_PER_NODE=1

######################################################################################################






########################## SETUP SOME VARIABLES AND PARAMETERS (AUTOMATIC)  ########################

# set the code BASE_DIR and the OUTPUT_BASE_DIR
if [ "$COMPUTER_TYPE" = "laptop" ] ; then
	BASE_DIR="/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes"
	OUTPUT_BASE_DIR="$BASE_DIR/results"
elif [ "$COMPUTER_TYPE" = "compute-lung" ] ; then
	BASE_DIR="/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes"
	OUTPUT_BASE_DIR="$BASE_DIR/results"
elif [ "$COMPUTER_TYPE" = "arcus-a" ] ; then
	BASE_DIR="/home/comp-respiratory-modelling/jmbewu/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes"
	OUTPUT_BASE_DIR="/home/comp-respiratory-modelling/jmbewu/data/results"
elif [ "$COMPUTER_TYPE" = "arcus-b" ] ; then
	BASE_DIR="/home/comp-respiratory-modelling/jmbewu/libmesh-git-b/libmesh/examples/dphil/coupled_navier_stokes"
	OUTPUT_BASE_DIR="/home/comp-respiratory-modelling/jmbewu/data/results"
else
	echo "ERROR: invalid COMPUTER_TYPE specified."
fi

# variables for script
THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
THIS_FILE="$0"



#############################################################################





############################ SIMULATION PARAMETERS ###########################

# simulation parameters
OUTPUT_DIR_RELATIVE="particle_deposition_august/new_output/particle_0d_unsteady_from_steady_real/"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$OUTPUT_DIR_RELATIVE"

MESH_BASE_DIR="$BASE_DIR"
MESH_BASE_DIR="/home/james/meshes-git/dphil_meshes"

CONSTANT_VARIABLES="--use-petsc --solver_variable_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_w 0 --solver_group_ns3d_p 1  --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_p 0 --solver_group_ns3d1d_w 0 --solver_group_ns3d1d_Q 1 --solver_group_ns3d1d_P 1 --solver_group_ns3d1d_0_u 0 --solver_group_ns3d1d_0_v 0 --solver_group_ns3d1d_0_w 0 --solver_group_ns3d1d_0_p 1 --solver_system_names -log_summary"

SEMI_CONSTANT_VARIABLES="-_read_1d_mesh true -num_1d_trees 1 -input_1d_file $MESH_BASE_DIR/APLE0036266/output_airways_full -radius_on_edge false -mesh_input_scaling_1d 1.e-3  -alveolated_1d_tree 0  -threed true -flow_mag_1d 150.0e-6 -viscosity 1.8e-5 -density 1.2 -stokes true -output_nondim false -match_1d_mesh_to_3d_mesh false -unsteady 0 -sim_type 1 -reynolds_number_calculation false -particle_deposition 5 -restart_folder $OUTPUT_BASE_DIR/particle_deposition_august/new_output/output_0d_real/ -unsteady_from_steady true -write_interval 0.001"

LIBMESH_OPTIONS="-input_file $THIS_DIR/navier.in -output_folder $OUTPUT_DIR/ $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES"

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

	echo "Submitted job $JOB_ID"
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
cp $THIS_DIR/$THIS_FILE $OUTPUT_DIR/

# output the git version in a file
GIT_VERSION=$(git describe)
echo "$GIT_VERSION" > $OUTPUT_DIR/git_version.dat

# output the name of the computer and the directory it was run from in a file
COMPUTER_NAME=$(hostname)
echo "$COMPUTER_NAME" > $OUTPUT_DIR/computer_name.dat
echo "$THIS_DIR" > $OUTPUT_DIR/run_directory.dat

##############################################################
