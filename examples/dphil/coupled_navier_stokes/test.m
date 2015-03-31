close all
A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/36266_new/0.01/deposition.dat'));

acinar_deposition_0p01 = A.data(:,10);
tb_deposition_0p01 = A.data(:,7);
total_deposition_0p01 = acinar_deposition_0p01 + tb_deposition_0p01;
subplot(1,4,1);
plot(acinar_deposition_0p01,'r-*')
hold on
plot(tb_deposition_0p01,'b-*')
plot(total_deposition_0p01,'g-*')
title('0.01\mu m')
xlabel('generation')

A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/36266_new/0.1/deposition.dat'));

acinar_deposition_0p1 = A.data(:,10);
tb_deposition_0p1 = A.data(:,7);
total_deposition_0p1 = acinar_deposition_0p1 + tb_deposition_0p1;
subplot(1,4,2);
plot(acinar_deposition_0p1,'r-*')
hold on
plot(tb_deposition_0p1,'b-*')
plot(total_deposition_0p1,'g-*')
title('0.1\mu m')
xlabel('generation')

A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/36266_new/1.0/deposition.dat'));

acinar_deposition_1 = A.data(:,10);
tb_deposition_1 = A.data(:,7);
total_deposition_1 = acinar_deposition_1 + tb_deposition_1;
subplot(1,4,3);
plot(acinar_deposition_1,'r-*')
hold on
plot(tb_deposition_1,'b-*')
plot(total_deposition_1,'g-*')
title('1\mu m')
xlabel('generation')

A = importdata(strcat('~/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/36266_new/10.0/deposition.dat'));

acinar_deposition_10 = A.data(:,10);
tb_deposition_10 = A.data(:,7);
total_deposition_10 = acinar_deposition_10 + tb_deposition_10;
subplot(1,4,4);
plot(acinar_deposition_10,'r-*')
hold on
plot(tb_deposition_10,'b-*')
plot(total_deposition_10,'g-*')
title('10\mu m')


figure
subplot(1,3,1);
total_deposition = [sum(total_deposition_0p01) sum(total_deposition_0p1) sum(total_deposition_1) sum(total_deposition_10)];
particle_size = [0.01 0.1 1.0 10.0];
semilogx(particle_size,total_deposition,'*-');
title('total dep - p-size')
subplot(1,3,2);
acinar_deposition = [sum(acinar_deposition_0p01) sum(acinar_deposition_0p1) sum(acinar_deposition_1) sum(acinar_deposition_10)];
particle_size = [0.01 0.1 1.0 10.0];
semilogx(particle_size,acinar_deposition,'*-');
title('acinar dep - p-size')
subplot(1,3,3);
tb_deposition = [sum(tb_deposition_0p01) sum(tb_deposition_0p1) sum(tb_deposition_1) sum(tb_deposition_10)];
particle_size = [0.01 0.1 1.0 10.0];
semilogx(particle_size,tb_deposition,'*-');
title('tb dep - p-size')