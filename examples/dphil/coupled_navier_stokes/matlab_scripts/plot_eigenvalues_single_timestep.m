% plot the eigenvalue spectrum of the preconditioned operator


A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/resistance_scaling/test6/eigenvalues.dat'));
num_time_steps = 100;
time_steps = A.data(:,1);
nonlinear_steps = A.data(:,2);
real_eigs = A.data(:,3:2:end);
imag_eigs = A.data(:,4:2:end);


max_nonlin = max(nonlinear_steps);
max_timestep = time_steps(end);

max_real = max(max(real_eigs));
min_real = min(min(real_eigs));
max_imag = max(max(imag_eigs));
min_imag = min(min(imag_eigs));

%close all
figure
plot(real_eigs,imag_eigs,'*');
axis([min_real max_real min_imag max_imag]);


