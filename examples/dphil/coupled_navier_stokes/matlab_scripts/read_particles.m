filename = '/home/james/dphil/analysis/particle_deposition/impaction/bifurcation/straight/1000_100/particles4000.dat';
A = importdata(filename);

particle_data = A.data;

num_particles = size(particle_data,1);

num_on_wall = 0;
num_on_surface = 0;

lost_on_wall = 0;

z_coord = particle_data(:,6);
exit_surface = particle_data(:,7);
exit_time = particle_data(:,8);
entrance_time = particle_data(:,9);

for i=1:num_particles
    if(exit_time(i) > 0)
        if(exit_surface(i) == -1)
            num_on_wall = num_on_wall + 1;
        %lost particles
        elseif(exit_surface(i) == -2)
            %if particle got lost above z=-3 then it's probably got lost
            %getting on the surface
            if(z_coord(i) > -0.03)
                num_on_wall = num_on_wall +1;
                lost_on_wall = lost_on_wall + 1;
            else
                num_on_surface = num_on_surface + 1;
            end
        %particles on exit surface
        elseif(exit_surface(i) > -1)
            num_on_surface = num_on_surface + 1;
        end
    end        
end

num_on_wall
lost_on_wall
num_on_surface