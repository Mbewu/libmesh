lambda = 68e-9;
particle_diameter = 0.01e-6;
cunningham_correction_factor = 1. + 2.*lambda/particle_diameter*(1.257 + 0.4*exp(-1.1*(particle_diameter/(2*lambda))))