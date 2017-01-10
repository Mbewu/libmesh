#!/bin/bash

# this is an example script that runs two simulations
# on a batch system it should definitely just submit them all, let's go with this assumption
# a script to run each simulation is created in that directory and then run

#########
# user should edit "some general user parameters" and of course the parameters for each libmesh simulation
#  similar to old scripts.





#############################################################################



########### SOME GENERAL USER PARAMETERS #############
COMPUTER_TYPE="arcus-b"
NUM_PROCS=1
QUEUE_TYPE="devel"
JOB_NAME="arc-job"
WALLTIME="00:10:00"
NUM_NODES=1
NUM_PROCS_PER_NODE=1


########################## SETUP SOME VARIABLES AND PARAMETERS (AUTOMATIC)  ########################

# set the code BASE_DIR and the OUTPUT_BASE_DIR
if [ "$COMPUTER_TYPE" = "laptop" ] ; then
	BASE_DIR="/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes"
	OUTPUT_BASE_DIR="$BASE_DIR/results"
elif [ "$COMPUTER_TYPE" = "compute-lung" ] ; then
	BASE_DIR="/users/jmbewu/coupled_navier_stokes"
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

OUTER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTER_FILE="$0"

#####################################################################################

########### SOME GENERAL LIBMESH/PETSC PARAMETERS FOR ALL SIMULATIONS ###############
#MESH_BASE_DIR="$BASE_DIR/meshes"
MESH_BASE_DIR="~/meshes-git/dphil_meshes"

CONSTANT_VARIABLES="--use-petsc --solver_variable_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_w 0 --solver_group_ns3d_p 1  --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_p 0 --solver_group_ns3d1d_w 0 --solver_group_ns3d1d_Q 1 --solver_group_ns3d1d_P 1 --solver_group_ns3d1d_0_u 0 --solver_group_ns3d1d_0_v 0 --solver_group_ns3d1d_0_w 0 --solver_group_ns3d1d_0_p 1 --solver_system_names -log_summary"

SEMI_CONSTANT_VARIABLES="-mesh_file $MESH_BASE_DIR/meshes_for_paper/multi_bifurcation_symmetric_unstructured_2.msh -num_generations_string 11 -inflow_bdy_id 0 -no_refinement 0 -output_no_refinement 0 -preconditioner_type_3d1d 10 -preconditioner_type_3d 9 -unsteady 1 -dt 5.0 -write_interval 5.0 -end_time 10.0 -period 1000.0 -stokes false -reynolds_number 500 -sim_type 0 -twod_oned_tree true -threed false -_read_1d_mesh false -0d_bifurcation_angle 0.392699082 -multiple_column_solve 1  -newton 0 -streamline_diffusion true -stokes true -fieldsplit true -no_refinement 0 -output_no_refinement 0 -linear_shape_functions true -stab true -preconditioner_type_3d1d 12 -preconditioner_type_3d 9 -reuse_preconditioner false"


########## FIRST SIMULATION #########
OUTPUT_DIR_RELATIVE="test_multiscript/test_1"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$OUTPUT_DIR_RELATIVE"
RUN_DIR="$OUTER_DIR/test_1"	# the directory 

LIBMESH_OPTIONS="-input_file $OUTER_DIR/navier.in -output_folder $OUTPUT_DIR/ -sim_type 5 -dt 5.0  $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES"

# generate the test_script file
$BASE_DIR/generate_test_script_file.sh "$COMPUTER_TYPE" "$NUM_PROCS" "$QUEUE_TYPE" "$JOB_NAME" "$WALLTIME" "$NUM_NODES" "$NUM_PROCS_PER_NODE" "$OUTPUT_DIR" "$LIBMESH_OPTIONS" "$RUN_DIR" "$BASE_DIR"
# make it runnable
chmod +x $RUN_DIR/test_script.sh
# run the test_script (which will generate a job_script and run it etc
cd $RUN_DIR
$RUN_DIR/test_script.sh
cd $OUTER_DIR

# done running


########## SECOND SIMULATION #########
OUTPUT_DIR_RELATIVE="test_multiscript/test_2"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$OUTPUT_DIR_RELATIVE"
RUN_DIR="$OUTER_DIR/test_2"	# the directory 

LIBMESH_OPTIONS="-input_file $OUTER_DIR/navier.in -output_folder $OUTPUT_DIR/ -sim_type 5 -dt 2.5  $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES"

# generate the test_script file
$BASE_DIR/generate_test_script_file.sh "$COMPUTER_TYPE" "$NUM_PROCS" "$QUEUE_TYPE" "$JOB_NAME" "$WALLTIME" "$NUM_NODES" "$NUM_PROCS_PER_NODE" "$OUTPUT_DIR" "$LIBMESH_OPTIONS" "$RUN_DIR" "$BASE_DIR"
# make it runnable
chmod +x $RUN_DIR/test_script.sh
# run the test_script (which will generate a job_script and run it etc
cd $RUN_DIR
$RUN_DIR/test_script.sh
cd $OUTER_DIR

# done running


################# OUTPUT SOME HOUSEKEEPING DATA #####################

# copy this script to the output folder
echo hi
cp $OUTER_DIR/$OUTER_FILE "$OUTPUT_BASE_DIR/test_multiscript"
echo bye

##############################################################

