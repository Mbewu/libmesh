% file that extracts data from the linear_iterations.dat files for the
% exact_preconditioner study.

% in the exact preconditioner:
% three preconditioners:        diagonal, schur_0d and schur_stokes
% two types bcs:                asymmetric or symmetric geometries
% three numbers of outlets:     2, 8 and 32 outlets
% three timestep sizes:         0.1,0.01,0.001


preconditioner_type = cellstr(['monolithic/ ';'partitioned/']);
bc_type = cellstr(['asym/'; 'sym/ ']);
num_outlets = cellstr(['outlets_2/ '; 'outlets_8/ '; 'outlets_32/']);
time_steps = cellstr(['dt_0.1/  '; 'dt_0.01/ '; 'dt_0.005/'; 'dt_0.001/']);
time_steps_double = [0.1;0.01;0.005;0.001];

base_folder = '~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/results_for_paper/results_for_paper_linear_iterations/monolithic_motivation/';

% diagonal, asym, 2, 0.1
    for prec=1:2

    file_id = [prec,1,3,1];


    filename = [base_folder,preconditioner_type{file_id(1)},bc_type{file_id(2)},num_outlets{file_id(3)},time_steps{file_id(4)},'linear_iterations.dat'];
    disp(['Checking: ',filename]);
    disp(['Preconditioner type ',preconditioner_type{file_id(1)}]);
    disp(['BC type ',bc_type{file_id(2)}]);
    disp(['Num outlets ',num_outlets{file_id(3)}]);
    disp(['Time step ',time_steps{file_id(4)}]);
    if(exist(filename, 'file') == 2)
        A = importdata(filename);

        total_linear_iterations = A.data(end,3);
        max_linear_iterations = A.data(end,4);
        total_nonlinear_iterations = size(A.data,1);
        total_time_steps = A.data(end,1);
        average_nonlinear_iterations = total_nonlinear_iterations/total_time_steps;
        average_linear_iterations_per_time_step = total_linear_iterations/total_time_steps;
        average_linear_iterations_per_nonlinear_iteration = total_linear_iterations/total_nonlinear_iterations;
        
        
        disp(['tot lin its: ',num2str(total_linear_iterations)]);
        disp(['max lin its: ',num2str(max_linear_iterations)]);
        disp(['tot nlin its: ',num2str(total_nonlinear_iterations)]);
        disp(['ave nlin its: ',num2str(average_nonlinear_iterations)]);
        disp(['ave lin its per timestep: ',num2str(average_linear_iterations_per_time_step)]);
        disp(['ave lin its per nlin it: ',num2str(average_linear_iterations_per_nonlinear_iteration)]);

        expected_time_steps = 1/time_steps_double(file_id(4));
        
        if(total_time_steps < expected_time_steps)
            if(A.data(end,7) - A.data(end-1,7) > 0)
                disp('Error: simulation did not complete because nonlinear divergence.');
            else
                disp('Error: simulation did not complete.');
            end
        else
            
            if(max(A.data(:,2)) > 10)
                disp('Error: simulation completes, but some timesteps took > 10 picard iterations.');
            end
        end
    else
        disp('Error: file does not exist.');    
    end
end
