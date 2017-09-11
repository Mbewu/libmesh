% function to calculate the yeh impaction probability
% St - Stokes number, possibly a vector over all generations
% branching_angle - branching angle (rad), possibly a vector over all generations
% 
function [p_imp_yeh] = yeh_impaction(St,branching_angle)

p_imp_yeh = zeros(1,length(St));
for i=1:length(St)
    if branching_angle(i)*St(i) < 1
        p_imp_yeh(i) = 1 - 2/pi.*acos(branching_angle(i)*St(i)) + 1./pi*sin(2.*acos(branching_angle(i)*St(i)));
    else
        p_imp_yeh(i) = 1;
    end
end

end