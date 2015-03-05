#!/bin/bash


#mpirun -np 1 ./example-opt -output_folder results/coupling_study/monolithic/0.001/ -dt 0.001 -sim_type  5
mpirun -np 1 ./example-opt -output_folder results/coupling_study/monolithic/0.02/ -dt 0.02 -sim_type  5
mpirun -np 1 ./example-opt -output_folder results/coupling_study/monolithic/0.1/ -dt 0.1 -sim_type  5
mpirun -np 1 ./example-opt -output_folder results/coupling_study/monolithic/0.2/ -dt 0.2 -sim_type  5

#mpirun -np 1 ./example-opt -output_folder results/coupling_study/explicit/0.001/ -dt 0.001 -sim_type  4
mpirun -np 1 ./example-opt -output_folder results/coupling_study/explicit/0.02/ -dt 0.02 -sim_type  4
mpirun -np 1 ./example-opt -output_folder results/coupling_study/explicit/0.1/ -dt 0.1 -sim_type  4
mpirun -np 1 ./example-opt -output_folder results/coupling_study/explicit/0.2/ -dt 0.2 -sim_type  4

#mpirun -np 1 ./example-opt -output_folder results/coupling_study/implicit/0.001/ -dt 0.001 -sim_type  3
mpirun -np 1 ./example-opt -output_folder results/coupling_study/implicit/0.02/ -dt 0.02 -sim_type  3
mpirun -np 1 ./example-opt -output_folder results/coupling_study/implicit/0.1/ -dt 0.1 -sim_type  3
mpirun -np 1 ./example-opt -output_folder results/coupling_study/implicit/0.2/ -dt 0.2 -sim_type  3




