#!/bin/bash

#mkdir results/coupling_study_new_new/
#mkdir results/coupling_study_new_new/reynolds_1000/
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.1/
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.1/monolithic
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.1/monolithic/ -dt 100.0 -sim_type 5 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.1/implicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.1/implicit/ -dt 100.0 -sim_type 3 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.1/explicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.1/explicit/ -dt 100.0 -sim_type 4 -reynolds_number 1000 -period 2000.0 -end_time 1000.0

#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.05/
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.05/monolithic
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.05/monolithic/ -dt 50.0 -sim_type 5 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.05/implicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.05/implicit/ -dt 50.0 -sim_type 3 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.05/explicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.05/explicit/ -dt 50.0 -sim_type 4 -reynolds_number 1000 -period 2000.0 -end_time 1000.0


#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.01/
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.01/monolithic
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.01/monolithic/ -dt 10.0 -sim_type 5 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.01/implicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.01/implicit/ -dt 10.0 -sim_type 3 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
#mkdir results/coupling_study_new_new/reynolds_1000/dt_0.01/explicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.01/explicit/ -dt 10.0 -sim_type 4 -reynolds_number 1000 -period 2000.0 -end_time 1000.0

mkdir results/coupling_study_new_new/reynolds_1000/dt_0.005/
mkdir results/coupling_study_new_new/reynolds_1000/dt_0.005/monolithic
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.005/monolithic/ -dt 5.0 -sim_type 5 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
mkdir results/coupling_study_new_new/reynolds_1000/dt_0.005/implicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.005/implicit/ -dt 5.0 -sim_type 3 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
mkdir results/coupling_study_new_new/reynolds_1000/dt_0.005/explicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new_new/reynolds_1000/dt_0.005/explicit/ -dt 5.0 -sim_type 4 -reynolds_number 1000 -period 2000.0 -end_time 1000.0
