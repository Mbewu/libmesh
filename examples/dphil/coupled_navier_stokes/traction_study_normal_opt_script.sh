#!/bin/bash
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/normal/25/ -velocity_scale 25 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/normal/15/ -velocity_scale 15 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/normal/5/ -velocity_scale 5 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/normal/1/ -velocity_scale 1 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/normal/0.1/ -velocity_scale 0.1 -optimisation_stabilised  1 -optimisation_scaling_factor 0

mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/25/ -optimisation_stabilised  1 -volume_optimisation  true -velocity_scale 25 
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/15/ -optimisation_stabilised  1 -volume_optimisation  true -velocity_scale 15
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/5/ -optimisation_stabilised  1 -volume_optimisation  true -velocity_scale 5
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/1/ -optimisation_stabilised  1 -volume_optimisation  true -velocity_scale 1
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/0.1/ -optimisation_stabilised  1 -volume_optimisation  true -velocity_scale 0.1


mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/25/ -optimisation_stabilised  1 -tangential_optimisation  true -velocity_scale 25
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/15/ -optimisation_stabilised  1 -tangential_optimisation  true -velocity_scale 15
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/5/ -optimisation_stabilised  1 -tangential_optimisation  true -velocity_scale 5
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/1/ -optimisation_stabilised  1 -tangential_optimisation  true -velocity_scale 1
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/0.1/ -optimisation_stabilised  1 -tangential_optimisation  true -velocity_scale 0.1



mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/25/ -neumann_stabilised true -velocity_scale 25 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/15/ -neumann_stabilised true -velocity_scale 15 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/5/ -neumann_stabilised true -velocity_scale 5 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/1/ -neumann_stabilised true -velocity_scale 1 -optimisation_stabilised  1 -optimisation_scaling_factor 0
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/0.1/ -neumann_stabilised true -velocity_scale 0.1 -optimisation_stabilised  1 -optimisation_scaling_factor 0
