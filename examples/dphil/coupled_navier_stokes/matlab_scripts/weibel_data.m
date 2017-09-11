% function that returns weibel symmetric tree and yeh data
function [diameter, branch_length, branching_angle, gravity_angle, generations] = weibel_data()

generations = 0:23;
diameter = [1.8 1.22 0.83 0.56 0.45 0.35 0.28 0.23 0.186 0.154 0.13 0.109 ...
    0.095 0.082 0.074 0.066 0.06 0.054 0.05 0.047 0.045 0.043 0.041 0.041]; %cm
branch_length = [12 4.76 1.9 0.76 1.27 1.07 0.9 0.76 0.64 0.54 0.46 0.39 0.33 0.27 ...
    0.23 0.2 0.165 0.141 0.117 0.099 0.083 0.07 0.059 0.05]; % cm

%yeh 1980 angles
branching_angle = pi/180*[0 33 34 22 20 18 19 22 28 22 33 34 37 39 39 51 45 45 45 45 45 45 45 45];
gravity_angle =   pi/180*[0 20 31 43 39 39 40 36 39 45 43 45 45 60 60 60 60 60 60 60 60 60 60 60];
   

end