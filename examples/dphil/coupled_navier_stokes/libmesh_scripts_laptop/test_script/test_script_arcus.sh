#!/bin/bash

# in this file we will run asymmetric simulations 2 outlets, dt=0.1,0.01,0.001, 8 outlets dt=0.1,0.01,0.001

# some variables
THIS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
THIS_FILE="test_script_1.sh"
BASE_DIR="/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes"


echo "The script is called $THIS_DIR/$THIS_FILE"

CONSTANT_VARIABLES="--use-petsc --solver_variable_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_w 0 --solver_group_ns3d_p 1  --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_p 0 --solver_group_ns3d1d_w 0 --solver_group_ns3d1d_Q 1 --solver_group_ns3d1d_P 1 --solver_group_ns3d1d_0_u 0 --solver_group_ns3d1d_0_v 0 --solver_group_ns3d1d_0_w 0 --solver_group_ns3d1d_0_p 1 --solver_system_names -log_summary -newton 0 -streamline_diffusion true -stokes true -fieldsplit true -no_refinement 0 -output_no_refinement 0 -linear_shape_functions true -stab true -preconditioner_type_3d1d 12 -preconditioner_type_3d 9 -reuse_preconditioner false"

SEMI_CONSTANT_VARIABLES="-reynolds_number 500 -period 1000.0 -end_time 500.0 -mesh_file $BASE_DIR/meshes/meshes_for_paper/multi_bifurcation_symmetric_unstructured_2.msh -num_generations_string 99 -unsteady 0"

# 2 dt 0.01
OUTPUT_FOLDER="$BASE_DIR/results/reproducibility_testing/test_1"
# command line options
LIBMESH_OPTIONS="-input_file $THIS_DIR/navier.in -output_folder $OUTPUT_FOLDER/ -sim_type 5 -dt 5.0  $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES"

# run the program
$BASE_DIR/example-opt -input_file $THIS_DIR/navier.in -output_folder $OUTPUT_FOLDER/ -sim_type 5 -dt 5.0  $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES 2>&1 | tee $OUTPUT_FOLDER/output.log
# copy this script to the output folder
cp $THIS_DIR/$THIS_FILE $OUTPUT_FOLDER/
# output the git version in a file
GIT_VERSION=$(git describe)
echo "$GIT_VERSION" > $OUTPUT_FOLDER/git_version.dat
# output the name of the computer in a file
COMPUTER_NAME=$(hostname)
echo "$COMPUTER_NAME" > $OUTPUT_FOLDER/computer_name.dat
echo "$THIS_DIR" > $OUTPUT_FOLDER/run_directory.dat

