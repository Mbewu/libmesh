
study_name='2d_lid_pcd_pre_p1p1';
Re=[10 50 100 200]% 1000];
size=[16 32 64 128];
%Re=[1];
%size=[16];
total_nonlin_gmres_iterations = zeros(length(Re),length(size));
stokes_gmres_iterations = zeros(length(Re),length(size));
max_gmres_iterations = zeros(length(Re),length(size));
nonlinear_iterations = zeros(length(Re),length(size));

for i=1:length(Re)
    for j=1:length(size)
        (strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/', ...
                study_name,'/Re',num2str(Re(i)),'/size',num2str(size(j)),'/linear_iterations.dat'))
        A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/', ...
                study_name,'/Re',num2str(Re(i)),'/size',num2str(size(j)),'/linear_iterations.dat'));
        total_nonlin_gmres_iterations(i,j) = A.data(end,3);
        stokes_gmres_iterations(i,j) = A.data(end,5);
        max_gmres_iterations(i,j) = A.data(end,4);
        nonlinear_iterations(i,j) = A.data(end,2);
    end
end

total_nonlin_gmres_iterations
stokes_gmres_iterations
max_gmres_iterations
average_nonlin_gmres_iterations = total_nonlin_gmres_iterations./(nonlinear_iterations-1)
nonlinear_iterations