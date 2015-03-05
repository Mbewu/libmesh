A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/1.0comp/norms.dat');

norms = zeros(size(A,1),9);
norms(:,1) = A(:,3);

A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/10.0comp/norms.dat');
norms(:,2) = A(:,3);
A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/100.0comp/norms.dat');
norms(:,3) = A(:,3);
A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/1000.0comp/norms.dat');
norms(:,4) = A(:,3);
A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/0.1comp/norms.dat');
norms(:,5) = A(:,3);
A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/0.01comp/norms.dat');
norms(:,6) = A(:,3);
A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/0.001comp/norms.dat');
norms(:,7) = A(:,3);
A = importdata('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/vol_opt_param_study/level_2/0.0001comp/norms.dat');
norms(:,8) = A(:,3);

time = 0:80;
time = time * 0.05; %zeros(size(A,1),1);



plot(time,norms)
legend('1.0','10.0','100.0','1000.0','0.1','0.01','0.001','0.0001');