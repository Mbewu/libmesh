% plot initial distribution of particles for a 2D simulation

filename = '/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/particle_deposition_august/test_particles/bifurcation_no_motion/1000_0_var/particles0000.dat';
A = importdata(filename);

particle_data = A.data;

num_particles = size(particle_data,1);

y_coord = particle_data(:,5);

hist(y_coord)
xlabel('position (cm)');
ylabel('num particles');
title('1000 particles, zero var');


