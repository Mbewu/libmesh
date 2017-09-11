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
OUTER_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTER_FILE="$0"


#############################################################################





############################ SIMULATION PARAMETERS ###########################


MESH_BASE_DIR="$BASE_DIR"
MESH_BASE_DIR="/home/james/meshes-git/dphil_meshes"

CONSTANT_VARIABLES="--use-petsc --solver_variable_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_w 0 --solver_group_ns3d_p 1  --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_p 0 --solver_group_ns3d1d_w 0 --solver_group_ns3d1d_Q 1 --solver_group_ns3d1d_P 1 --solver_group_ns3d1d_0_u 0 --solver_group_ns3d1d_0_v 0 --solver_group_ns3d1d_0_w 0 --solver_group_ns3d1d_0_p 1 --solver_system_names -log_summary -preconditioner_type_3d 0 -preconditioner_type_3d1d 0"

SEMI_CONSTANT_VARIABLES="-geometry_type 4 -mesh_file $MESH_BASE_DIR/meshes_for_deposition/multi_bifurcation_symmetric_diff_wall_bdy_2.msh -mesh_input_scaling_3d 0.02 -inflow_bdy_id 1011 -gmsh_diff_wall_bdy_id true -_read_1d_mesh false -num_generations_string 0 -num_generations 8 -twod_oned_tree true -threed false -stokes false -unsteady 1 -sim_type 3 -match_1d_to_3d_mesh true -reynolds_number_calculation false -prescribed_flow 0 -period 4.0 -end_time 2.0 -dt 0.01 -write_interval 0.01 -nonlinear_tolerance 1e-6 -density 1.2 -viscosity 1.8e-5 -newton 1 -bifurcation_start_1d true -half_initial_length false -use_centreline_data true -input_1d_file  $MESH_BASE_DIR/meshes_for_deposition/multi_bifurcation_symmetric_diff_wall_bdy_2 -radius_on_edge false -particle_deposition 6 -unsteady_from_steady false -restart_folder $OUTPUT_BASE_DIR/particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_30/fluid_sim/ -velocity_mag_3d 3.183 -particle_end_time 2.0 -particle_deposition_rate 10000 -particle_dt 0.00003333333333333333333333333333333"


###########################################################################








####################### FIRST SIMULATION ##########################

# simulation parameters
OUTPUT_DIR_RELATIVE="particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_30/particle_size_2/"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$OUTPUT_DIR_RELATIVE"
RUN_DIR="$OUTER_DIR/particle_size_2/"	# the directory 

LIBMESH_OPTIONS="-input_file $OUTER_DIR/navier.in -particle_deposition_input_file $OUTER_DIR/particle_deposition.in -particle_diameter 2.e-6 -output_folder $OUTPUT_DIR/ $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES"



# generate the test_script file
$BASE_DIR/generate_test_script_file.sh "$COMPUTER_TYPE" "$NUM_PROCS" "$QUEUE_TYPE" "$JOB_NAME" "$WALLTIME" "$NUM_NODES" "$NUM_PROCS_PER_NODE" "$OUTPUT_DIR" "$LIBMESH_OPTIONS" "$RUN_DIR" "$BASE_DIR"
# make it runnable
chmod +x $RUN_DIR/test_script.sh
# run the test_script (which will generate a job_script and run it etc
cd $RUN_DIR
$RUN_DIR/test_script.sh
cd $OUTER_DIR

# simulation finished running

###########################################################################


####################### SECOND SIMULATION ##########################

# simulation parameters
OUTPUT_DIR_RELATIVE="particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_30/particle_size_5/"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$OUTPUT_DIR_RELATIVE"
RUN_DIR="$OUTER_DIR/particle_size_5/"	# the directory 

LIBMESH_OPTIONS="-input_file $OUTER_DIR/navier.in -particle_deposition_input_file $OUTER_DIR/particle_deposition.in -particle_diameter 5.e-6 -output_folder $OUTPUT_DIR/ $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES"



# generate the test_script file
$BASE_DIR/generate_test_script_file.sh "$COMPUTER_TYPE" "$NUM_PROCS" "$QUEUE_TYPE" "$JOB_NAME" "$WALLTIME" "$NUM_NODES" "$NUM_PROCS_PER_NODE" "$OUTPUT_DIR" "$LIBMESH_OPTIONS" "$RUN_DIR" "$BASE_DIR"
# make it runnable
chmod +x $RUN_DIR/test_script.sh
# run the test_script (which will generate a job_script and run it etc
cd $RUN_DIR
$RUN_DIR/test_script.sh
cd $OUTER_DIR

# simulation finished running

###########################################################################


####################### THIRD SIMULATION ##########################

# simulation parameters
OUTPUT_DIR_RELATIVE="particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_30/particle_size_10/"
OUTPUT_DIR="$OUTPUT_BASE_DIR/$OUTPUT_DIR_RELATIVE"
RUN_DIR="$OUTER_DIR/particle_size_10/"	# the directory 

LIBMESH_OPTIONS="-input_file $OUTER_DIR/navier.in -particle_deposition_input_file $OUTER_DIR/particle_deposition.in -particle_diameter 10.e-6 -output_folder $OUTPUT_DIR/ $CONSTANT_VARIABLES $SEMI_CONSTANT_VARIABLES"



# generate the test_script file
$BASE_DIR/generate_test_script_file.sh "$COMPUTER_TYPE" "$NUM_PROCS" "$QUEUE_TYPE" "$JOB_NAME" "$WALLTIME" "$NUM_NODES" "$NUM_PROCS_PER_NODE" "$OUTPUT_DIR" "$LIBMESH_OPTIONS" "$RUN_DIR" "$BASE_DIR"
# make it runnable
chmod +x $RUN_DIR/test_script.sh
# run the test_script (which will generate a job_script and run it etc
cd $RUN_DIR
$RUN_DIR/test_script.sh
cd $OUTER_DIR

# simulation finished running

###########################################################################



################# OUTPUT SOME HOUSEKEEPING DATA #####################

# copy this script to the output folder
echo hi
cp $OUTER_DIR/$OUTER_FILE "$OUTPUT_BASE_DIR/test_multiscript"
echo bye

##############################################################

