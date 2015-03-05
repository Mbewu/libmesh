#!/bin/bash
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/25/ -results_folder_1 results/stab_bc_study/normal/25/ -results_folder_2 results/stab_bc_study/traction/25/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/15/ -results_folder_1 results/stab_bc_study/normal/15/ -results_folder_2 results/stab_bc_study/traction/15/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/5/ -results_folder_1 results/stab_bc_study/normal/5/ -results_folder_2 results/stab_bc_study/traction/5/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/1/ -results_folder_1 results/stab_bc_study/normal/1/ -results_folder_2 results/stab_bc_study/traction/1/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/traction/0.1/ -results_folder_1 results/stab_bc_study/normal/0.1/ -results_folder_2 results/stab_bc_study/traction/0.1/ -compare_results true


mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/25/ -results_folder_1 results/stab_bc_study/normal_opt/25/ -results_folder_2 results/stab_bc_study/vol_opt/25/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/15/ -results_folder_1 results/stab_bc_study/normal_opt/15/ -results_folder_2 results/stab_bc_study/vol_opt/15/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/5/ -results_folder_1 results/stab_bc_study/normal_opt/5/ -results_folder_2 results/stab_bc_study/vol_opt/5/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/1/ -results_folder_1 results/stab_bc_study/normal_opt/1/ -results_folder_2 results/stab_bc_study/vol_opt/1/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/vol_opt/0.1/ -results_folder_1 results/stab_bc_study/normal_opt/0.1/ -results_folder_2 results/stab_bc_study/vol_opt/0.1/ -compare_results true

mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/25/ -results_folder_1 results/stab_bc_study/normal_opt/25/ -results_folder_2 results/stab_bc_study/surf_opt/25/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/15/ -results_folder_1 results/stab_bc_study/normal_opt/15/ -results_folder_2 results/stab_bc_study/surf_opt/15/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/5/ -results_folder_1 results/stab_bc_study/normal_opt/5/ -results_folder_2 results/stab_bc_study/surf_opt/5/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/1/ -results_folder_1 results/stab_bc_study/normal_opt/1/ -results_folder_2 results/stab_bc_study/surf_opt/1/ -compare_results true
mpirun -np 1 ./example-opt -output_folder results/stab_bc_study/surf_opt/0.1/ -results_folder_1 results/stab_bc_study/normal_opt/0.1/ -results_folder_2 results/stab_bc_study/surf_opt/0.1/ -compare_results true
