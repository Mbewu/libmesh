#!/bin/bash

study_name=3d_bifurcpipeflow_pressure_pre_stokes_p2p1
fieldsplit_options="--solver_variable_names --solver_group_ns3d_u 0 --solver_group_ns3d_v 0 --solver_group_ns3d_p 1 --solver_group_ns3d_w 0 --solver_system_names ns3d"

# array of problem sizes
size=(0)
reynolds=(1 100 1000)

mkdir results/$study_name

for i in "${reynolds[@]}"; do
	echo itemRe: $i

	mkdir results/$study_name/Re$i/
	for j in "${size[@]}"; do

		mkdir results/$study_name/Re$i/size$j/
		mpirun -np 1 ./example-opt $fieldsplit_options -output_folder results/$study_name/Re$i/size$j/ -reynolds_number $i -no_refinement $j
		echo itemSize: $j
	done
done

#mkdir results/coupling_study_new/dt_0.1/
#mkdir results/coupling_study_new/dt_0.1/monolithic
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.1/monolithic/ -dt 0.1 -sim_type 5
#mkdir results/coupling_study_new/dt_0.1/implicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.1/implicit/ -dt 0.1 -sim_type 3
#mkdir results/coupling_study_new/dt_0.1/explicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.1/explicit/ -dt 0.1 -sim_type 4

#mkdir results/coupling_study_new/dt_0.01/
#mkdir results/coupling_study_new/dt_0.01/monolithic
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.01/monolithic/ -dt 0.01 -sim_type 5
#mkdir results/coupling_study_new/dt_0.01/implicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.01/implicit/ -dt 0.01 -sim_type 3
#mkdir results/coupling_study_new/dt_0.01/explicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.01/explicit/ -dt 0.01 -sim_type 4


#mkdir results/coupling_study_new/dt_0.002/
#mkdir results/coupling_study_new/dt_0.002/monolithic
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.002/monolithic/ -dt 0.002 -sim_type 5
#mkdir results/coupling_study_new/dt_0.002/implicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.002/implicit/ -dt 0.002 -sim_type 3
#mkdir results/coupling_study_new/dt_0.002/explicit
#mpirun -np 1 ./example-opt -output_folder results/coupling_study_new/dt_0.002/explicit/ -dt 0.002 -sim_type 4
