#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=picard_navier_stokes_read_mesh_1

# so that petsc is always used...
options="--use-petsc"
run_example "$example_name" "$options"
