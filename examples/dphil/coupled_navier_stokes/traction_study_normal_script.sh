#!/bin/bash
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/25/ -velocity_scale 25 -optimisation_stabilised  0 -optimisation_scaling_factor 0 -neumann_stabilised true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/15/ -velocity_scale 15 -optimisation_stabilised  0 -optimisation_scaling_factor 0 -neumann_stabilised true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/5/ -velocity_scale 5 -optimisation_stabilised  0 -optimisation_scaling_factor 0 -neumann_stabilised true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/1/ -velocity_scale 1 -optimisation_stabilised  0 -optimisation_scaling_factor 0 -neumann_stabilised true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/0.1/ -velocity_scale 0.1 -optimisation_stabilised  0 -optimisation_scaling_factor 0 -neumann_stabilised true

