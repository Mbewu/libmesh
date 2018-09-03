% function that returns weibel symmetric tree and yeh data
function [diameter, branch_length, branching_angle, gravity_angle, generations] = james_symmetric_data(max_generation,bif_angle,simple_gravity_angle)

gen_ratio = 0.8;
bif_angle = bif_angle;
d_start = 2.0;
l_start = 6.5;

generations = 0:max_generation;
diameter = d_start*gen_ratio.^generations;          % cm
branch_length = l_start*gen_ratio.^generations;     % cm

%yeh 1980 angles
branching_angle = pi/180.* bif_angle*ones(1,length(generations));

% need a way to calculate the average angle to gravity of each pipe
% can't think of anything except doing it separately for each generation
% TURNS OUT IT IS 45 for 45, BUT NOT FOR OTHER BIFURCATION ANGLES

if(simple_gravity_angle)
    gravity_angle = pi/180.* bif_angle*ones(1,length(generations));    
else
    gravity_angle = zeros(1,length(generations));

    gravity_angle(1) = 0;   % trachea
    gravity_angles_previous = 0;    % the trachea
    count = 1;  % the position we put stuff in
    for i=generations(2:end)
        count = count + 1;
        generation = i
        gravity_angles = zeros(2.^i,1);

        %loop over the gravity_angles_previous and populate gravity_angles
        for j=1:length(gravity_angles_previous)
            gravity_angles(2*j - 1) = gravity_angles_previous(j)+1;
            gravity_angles(2*j) = gravity_angles_previous(j)-1;
        end

        gravity_angles = gravity_angles
        % save the integer values before calculating the actual angles and
        % periodicity
        gravity_angles_previous = gravity_angles;

        % calculate actual angle
        gravity_angles = gravity_angles * bif_angle;

        % remove periodicity, the function will correctly handle negative
        % angles and give them positive values
        gravity_angles = mod(gravity_angles,360);

        % to get the angle to gravity, mod 180 
        % and then is angle > 90, angle = 180 - angle
        gravity_angles = mod(gravity_angles,180);
        gravity_angles(gravity_angles>90) = 180 - gravity_angles(gravity_angles>90)

        % average gravity angle
        gravity_angle(count) = mean(gravity_angles)
    end

    gravity_angle =   pi/180*gravity_angle
end

gravity_angle(1) = 0.;


end