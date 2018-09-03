% plot the eigenvalue spectrum of the preconditioned operator

exact = true;
approx = true;

step = 1;   % step to read

max_real = 0;
min_real = 0;
max_imag = 0;
min_imag = 0;

% the original source
%A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/resistance_scaling/test6/eigenvalues.dat'));
if(exact)
    A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/monolithic_paper_revisions/eigenvalues/mono_diag/R10/eigenvalues.dat'));
    num_time_steps = 100;
    time_steps = A.data(step,1);
    nonlinear_steps = A.data(step,2);
    real_eigs = A.data(step,3:2:end);
    imag_eigs = A.data(step,4:2:end);
    eigs = real_eigs + imag_eigs*i;



    max_nonlin = max(nonlinear_steps);
    max_timestep = time_steps(end);

    max_real = max(max(real_eigs));
    min_real = min(min(real_eigs));
    max_imag = max(max(imag_eigs));
    min_imag = min(min(imag_eigs));
end

if(approx)
    % there is a problem in that the data is not rectangular
    fid = fopen('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/monolithic_paper_revisions/eigenvalues/mono_diag/R10/eigenvalues_approx.dat');
    line_ex = fgetl(fid);
    line_ex = fgetl(fid); % skip first two lines
    % want to skip to line step
    line = 1;
    while line < step
        line_ex = fgetl(fid);
        line = line + 1;
    end
    %now we read the next line
    line_ex = fgetl(fid);
    numbers = str2num(line_ex);
    %chars = strsplit(line_ex);
    approx_real_eigs = numbers(3:2:end);
    approx_imag_eigs = numbers(4:2:end);
    %A_approx = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/monolithic_paper_revisions/eigenvalues/ns_no_pc/eigenvalues_approx.dat'));
    %approx_real_eigs = A_approx.data(step,3:2:end);
    %approx_imag_eigs = A_approx.data(step,4:2:end);
    approx_eigs = approx_real_eigs + approx_imag_eigs*i;

    max_approx_real = max(max(approx_real_eigs));
    min_approx_real = min(min(approx_real_eigs));
    max_approx_imag = max(max(approx_imag_eigs));
    min_approx_imag = min(min(approx_imag_eigs));
    
    if(exact)
        max_real = max(max_real,max_approx_real);
        min_real = max(min_real,min_approx_real);
        max_imag = max(max_imag,max_approx_imag);
        min_imag = max(min_imag,min_approx_imag);
    else
        max_real = max_approx_real;
        min_real = min_approx_real;
        max_imag = max_approx_imag;
        min_imag = min_approx_imag;
    end
end


limit_axes =false;
% always shift the limits to give space
max_real = 5;
min_real = -5;
max_imag = 1;
min_imag = -1;

%figure
%plot(real_eigs,imag_eigs,'*b','DisplayName','exact eigs');
%axis([min_real max_real min_imag max_imag]);

%close all
figure
subplot(1,3,1)
if(exact)
    plot(real_eigs,imag_eigs,'*b','DisplayName','exact eigs');
    if(limit_axes)
        axis([min_real max_real min_imag max_imag]);
    end
end
title('exact');

subplot(1,3,2)
if(approx)
    plot(approx_real_eigs,approx_imag_eigs,'*r','DisplayName','approx eigs');
    if(limit_axes)
        axis([min_real max_real min_imag max_imag]);
    end
end
title('approx');

subplot(1,3,3);
if(exact)
    plot(real_eigs,imag_eigs,'*b','DisplayName','exact eigs');
    if(limit_axes)
        axis([min_real max_real min_imag max_imag]);
    end
end
if(approx)
    hold on
    plot(approx_real_eigs,approx_imag_eigs,'*r','DisplayName','approx eigs');
    if(limit_axes)
        axis([min_real max_real min_imag max_imag]);
    end
end
title('combined');

if(exact)
    if(approx)
        legend('exact','approx');
    else
        legend('exact');
    end
elseif(approx)
	legend('approx');
end
    
    
%axis([min_real max_real min_imag max_imag]);


