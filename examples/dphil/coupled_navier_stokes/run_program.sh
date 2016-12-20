#!/bin/sh

############################################################
## SCRIPT FOR RUNNING LIBMESH PROGRAM ON NON-BATCH COMPUTER
## - parameters
## - $1 - BASE_DIR -> directory executable is in
## - $2 - LIBMESH_OPTIONS -> options passed to libmesh
## - $3 - OUTPUT_DIR -> output folder
## - $4 - NUM_PROCS -> number of procs
############################################################

BASE_DIR=$1
LIBMESH_OPTIONS=$2
OUTPUT_DIR=$3
NUM_PROCS=$4

# sequential
#$BASE_DIR/example-opt $LIBMESH_OPTIONS | tee $OUTPUT_DIR/output.log

# parallel (bash sucks with spaces...)
RUN_COMMAND="mpirun -np "$NUM_PROCS" "$BASE_DIR"/example-opt "$LIBMESH_OPTIONS""
$RUN_COMMAND | tee "$OUTPUT_DIR"/output.log

