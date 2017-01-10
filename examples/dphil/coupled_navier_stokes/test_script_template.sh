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

COMPUTER_TYPE="%COMPUTER_TYPE%"      #"laptop" is laptop, "compute-lung" is compute-lung, "arcus-a" is arcus-a and "arcus-b" is arcus-b

# non-batch specific variables
NUM_PROCS=%NUM_PROCS%

# batch specific variables
QUEUE_TYPE="%QUEUE_TYPE%"           #"devel" is test and "" is just normal
JOB_NAME="%JOB_NAME%"
WALLTIME=%WALLTIME%
NUM_NODES=%NUM_NODES%
NUM_PROCS_PER_NODE=%NUM_PROCS_PER_NODE%

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

OUTPUT_DIR="%OUTPUT_DIR%"
LIBMESH_OPTIONS="%LIBMESH_OPTIONS%"

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

