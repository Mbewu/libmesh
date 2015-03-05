#!/bin/bash

mkdir results/coupling_study_new/
mkdir results/coupling_study_new/dt_0.1/
mkdir results/coupling_study_new/dt_0.1/monolithic
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.1/monolithic/ -dt 0.1 -sim_type 5
mkdir results/coupling_study_new/dt_0.1/implicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.1/implicit/ -dt 0.1 -sim_type 3
mkdir results/coupling_study_new/dt_0.1/explicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.1/explicit/ -dt 0.1 -sim_type 4

mkdir results/coupling_study_new/dt_0.01/
mkdir results/coupling_study_new/dt_0.01/monolithic
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.01/monolithic/ -dt 0.01 -sim_type 5
mkdir results/coupling_study_new/dt_0.01/implicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.01/implicit/ -dt 0.01 -sim_type 3
mkdir results/coupling_study_new/dt_0.01/explicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.01/explicit/ -dt 0.01 -sim_type 4


mkdir results/coupling_study_new/dt_0.002/
mkdir results/coupling_study_new/dt_0.002/monolithic
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.002/monolithic/ -dt 0.002 -sim_type 5
mkdir results/coupling_study_new/dt_0.002/implicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.002/implicit/ -dt 0.002 -sim_type 3
mkdir results/coupling_study_new/dt_0.002/explicit
mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.002/explicit/ -dt 0.002 -sim_type 4
