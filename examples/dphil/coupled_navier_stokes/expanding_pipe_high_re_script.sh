#!/bin/bash
#mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/normal/
#mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/stabilised_0.2/ -neumann_stabilised 1 -backflow_stab_param 0.2
#mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/stabilised_1.0/ -neumann_stabilised 1 -backflow_stab_param 1.0
mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/optimised_test_1.0/ -optimisation_stabilised 1 -optimisation_scaling_factor 1.0 -tangential_optimisation 0 -volume_optimisation 1
mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/optimised_tangent/ -optimisation_stabilised 1 -optimisation_scaling_factor 1.0 -tangential_optimisation 1 -volume_optimisation 0
#mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/optimised_volume_0.0/ -optimisation_stabilised 1 -optimisation_scaling_factor 0.0 -tangential_optimisation 0 -volume_optimisation 1
#mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/optimised_volume_0.1/ -optimisation_stabilised 1 -optimisation_scaling_factor 0.1 -tangential_optimisation 0 -volume_optimisation 1
#mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/optimised_volume_0.01/ -optimisation_stabilised 1 -optimisation_scaling_factor 0.01 -tangential_optimisation 0 -volume_optimisation 1
#mpirun -np 1 ./example-opt -output_folder results/expanding_pipe_study_high_re/optimised_volume_10.0/ -optimisation_stabilised 1 -optimisation_scaling_factor 10.0 -tangential_optimisation 0 -volume_optimisation 1
