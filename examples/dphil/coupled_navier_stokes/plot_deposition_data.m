A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/deposition.dat'));
A

subplot(1,2,1)
plot(A.data(:,1),A.data(:,2))
title('generation plot')

subplot(1,2,2)
plot(A.data(:,3),A.data(:,4))
title('order plot')