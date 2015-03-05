#!/bin/bash

mkdir results/pressure_flow_study/poiseuille
mkdir results/pressure_flow_study/poiseuille/2.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/2.0/ -flow_mag_1d 2.0e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/1.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/1.75/ -flow_mag_1d 1.75e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/1.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/1.5/ -flow_mag_1d 1.5e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/1.25/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/1.25/ -flow_mag_1d 1.25e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/1.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/1.0/ -flow_mag_1d 1.0e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/0.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.75/ -flow_mag_1d 0.75e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/0.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.5/ -flow_mag_1d 0.5e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/0.4/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.4/ -flow_mag_1d 0.4e-3 -resistance_type_1d 
mkdir results/pressure_flow_study/poiseuille/0.3/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.3/ -flow_mag_1d 0.3e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/0.2/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.2/ -flow_mag_1d 0.2e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/0.1/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.1/ -flow_mag_1d 0.1e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/0.05/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.05/ -flow_mag_1d 0.05e-3 -resistance_type_1d 0
mkdir results/pressure_flow_study/poiseuille/0.01/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/poiseuille/0.01/ -flow_mag_1d 0.01e-3 -resistance_type_1d 0


mkdir results/pressure_flow_study/pedley
mkdir results/pressure_flow_study/pedley/2.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/2.0/ -flow_mag_1d 2.0e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/1.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/1.75/ -flow_mag_1d 1.75e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/1.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/1.5/ -flow_mag_1d 1.5e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/1.25/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/1.25/ -flow_mag_1d 1.25e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/1.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/1.0/ -flow_mag_1d 1.0e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.75/ -flow_mag_1d 0.75e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.5/ -flow_mag_1d 0.5e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.4/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.4/ -flow_mag_1d 0.4e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.3/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.3/ -flow_mag_1d 0.3e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.2/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.2/ -flow_mag_1d 0.2e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.1/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.1/ -flow_mag_1d 0.1e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.05/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.05/ -flow_mag_1d 0.05e-3 -resistance_type_1d 1
mkdir results/pressure_flow_study/pedley/0.01/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/pedley/0.01/ -flow_mag_1d 0.01e-3 -resistance_type_1d 1


mkdir results/pressure_flow_study/ertbruggen
mkdir results/pressure_flow_study/ertbruggen/2.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/2.0/ -flow_mag_1d 2.0e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/1.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/1.75/ -flow_mag_1d 1.75e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/1.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/1.5/ -flow_mag_1d 1.5e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/1.25/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/1.25/ -flow_mag_1d 1.25e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/1.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/1.0/ -flow_mag_1d 1.0e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.75/ -flow_mag_1d 0.75e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.5/ -flow_mag_1d 0.5e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.4/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.4/ -flow_mag_1d 0.4e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.3/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.3/ -flow_mag_1d 0.3e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.2/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.2/ -flow_mag_1d 0.2e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.1/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.1/ -flow_mag_1d 0.1e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.05/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.05/ -flow_mag_1d 0.05e-3 -resistance_type_1d 2
mkdir results/pressure_flow_study/ertbruggen/0.01/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/ertbruggen/0.01/ -flow_mag_1d 0.01e-3 -resistance_type_1d 2


mkdir results/pressure_flow_study/reynolds
mkdir results/pressure_flow_study/reynolds/2.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/2.0/ -flow_mag_1d 2.0e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/1.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/1.75/ -flow_mag_1d 1.75e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/1.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/1.5/ -flow_mag_1d 1.5e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/1.25/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/1.25/ -flow_mag_1d 1.25e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/1.0/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/1.0/ -flow_mag_1d 1.0e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.75/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.75/ -flow_mag_1d 0.75e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.5/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.5/ -flow_mag_1d 0.5e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.4/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.4/ -flow_mag_1d 0.4e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.3/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.3/ -flow_mag_1d 0.3e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.2/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.2/ -flow_mag_1d 0.2e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.1/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.1/ -flow_mag_1d 0.1e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.05/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.05/ -flow_mag_1d 0.05e-3 -resistance_type_1d 3
mkdir results/pressure_flow_study/reynolds/0.01/
mpirun -np 1 ./example-opt -output_folder results/pressure_flow_study/reynolds/0.01/ -flow_mag_1d 0.01e-3 -resistance_type_1d 3

