#!/bin/bash

#set -x

source $LIBMESH_DIR/examples/run_common.sh

example_name=picard_navier_stokes_read_mesh_1

# so that petsc is always used...
options="--use-petsc --solver_variable_names --solver_system_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_w 0 --solver_group_ns3d_p 1 --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_w 0 --solver_group_ns3d1d_p 0   --solver_group_ns3d1d_P 1  --solver_group_ns3d1d_Q 1  --solver_group_ns3d1d_0_u 0 --solver_group_ns3d1d_0_v 0  --solver_group_ns3d1d_0_w 0 --solver_group_ns3d1d_0_p 1 -log_summary -ns3d_pc_mg_log"
#options="--use-petsc --solver_variable_names --solver_group_ns3d1d_u 0 --solver_group_ns3d1d_v 0 --solver_group_ns3d1d_p 0 --solver_group_ns3d1d_w 0  --solver_group_ns3d1d_P 1  --solver_group_ns3d1d_Q 1 --solver_system_names ns3d1d -log_summary"
run_example "$example_name" "$options"
