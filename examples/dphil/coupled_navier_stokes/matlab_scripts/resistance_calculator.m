% calculate the resistance of a bifurcation tree

base_R = 1.;    % resistance of first pipe
scaling_factor = 0.85;  % 
num_bif = 2;    % number of bifurcation points

scale = 0.;
for i=0:num_bif
    scale = scale + 1./(scaling_factor^(3*i) * 2^i);
end

num_bif
total_resistance = base_R * scale