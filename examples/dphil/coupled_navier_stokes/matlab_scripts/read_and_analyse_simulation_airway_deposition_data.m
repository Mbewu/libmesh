% read in airway data and plot it

% read in 0d airway deposition data from last time step
folder = '/home/james/libmesh-git/libmesh/examples/dphil/coupled_navier_stokes/results/particle_deposition_august/parameter_sweep/first_try/geom_2/ifr_30/test/';
file = 'airway_deposition_data_3d_6000.dat';
full_path = [folder file];
airway_0d_deposition = importdata(full_path);
airway_0d_deposition_data = airway_0d_deposition.data;
