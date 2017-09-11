% diffusion coefficient (or particle diffusivity)
% stokes einstein equations

% D = k*T*C_s/(3*pi*mu*D_p)
k =  1.38064852e-23; % 1.38064852 × 10-23 boltzmann
T = 293.0;   % 20 C
mu = 1.8e-5; %air viscosity
lambda = 68.e-9; % mean free path
D_p = 1.0*1e-6;  % diameter
C_s = 1. + 2.*lambda./D_p.*(1.257 + 0.4*exp(-1.1*(D_p./(2*lambda)))) % cuningham

D = k*T*C_s/(3*pi*mu*D_p)
