#!/bin/sh

##############################################
## Script for creating PBS files on arcus-a from a template
## - parameters:
## - JOB_NAME - name of job (text)
## - WALLTIME - max walltime (00:00:00)
## - NUM_NODES - number of nodes
## - NUM_PROCS_PER_NODE - number of procs per node
## - LIBMESH_OPTIONS - command line options for libmesh program
## - OUTPUT_DIR - direrctory to output PBS script to
## - BASE_DIR - directory where application is
##############################################

OUTPUT_DIR="$6"
BASE_FOLDER="$7"

sed -e "s;%JOB_NAME%;$1;g" -e "s;%WALLTIME%;$2;g" -e "s;%NUM_NODES%;$3;g" -e "s;%NUM_PROCS_PER_NODE%;$4;g" -e "s;%LIBMESH_OPTIONS%;$5;g" -e "s;%BASE_DIR%;$7;g" $BASE_FOLDER/slurm_template.sh > $OUTPUT_DIR/job_script.sh
