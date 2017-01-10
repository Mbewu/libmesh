#!/bin/sh

##############################################
## Script for creating PBS files on arcus-a from a template
## - parameters:
## - COMPUTER_TYPE,1 - the type of computer you are running on (laptop, compute-lung, arcus-a or arcus-b)
## - NUM_PROCS,2 - number of procs, for non-batch computer
## - QUEUE_TYPE,3 - type of queue, devel for test queue, empty string for normal queue
## - JOB_NAME,4 - name of the job
## - WALLTIME,5 - max wall time in format "100:00:00"
## - NUM_NODES,6 - number of nodes on a batch system
## - NUM_PROCS_PER_NODE,7 - number of procs per node on a batch system
## - OUTPUT_DIR,8 - the directory for the results to go in
## - LIBMESH_OPTIONS,9 - the options to put to libmesh
##############################################

RUN_DIR="${10}"
BASE_FOLDER="${11}"


sed -e "s;%COMPUTER_TYPE%;$1;g" -e "s;%NUM_PROCS%;$2;g" -e "s;%QUEUE_TYPE%;$3;g" -e "s;%JOB_NAME%;$4;g" -e "s;%WALLTIME%;$5;g" -e "s;%NUM_NODES%;$6;g" -e "s;%NUM_PROCS_PER_NODE%;$7;g" -e "s;%OUTPUT_DIR%;$8;g" -e "s;%LIBMESH_OPTIONS%;$9;g" $BASE_FOLDER/test_script_template.sh > $RUN_DIR/test_script.sh
