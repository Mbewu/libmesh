% function to calculate the wang sedimentation probability
% St - Stokes number, possibly a vector over all generations
% gravity_angle - gravity angle, possibly a vector over all generations
% branch_length
% branch_velocity
function [p_sed_wang] = wang_sedimentation(St,gravity_angle,branch_length,branch_velocity)

g = 9.81;
p_sed_wang = zeros(1,length(St));
p_sed_wang = 1 - exp(-4*g*branch_length.*cos(gravity_angle).*St./(pi*branch_velocity.^2));

end



