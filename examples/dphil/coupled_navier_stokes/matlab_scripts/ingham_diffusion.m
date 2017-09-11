% function to calculate the yeh impaction probability
% diffusivity - particle diffusivity
% L - tube length
% Q - flow rate
function [p_dif_ingham] = ingham_diffusion(length,flow_rate,diffusivity)

diffusivity
p_dif_ingham = zeros(1,size(flow_rate,2));
%assuming flow rate > 0
for i=1:size(flow_rate,2)
    Delta = pi*diffusivity*length(i)/4/flow_rate(i);
    p_dif_ingham(i) = 1 - 0.819*exp(-14.63*Delta) - 0.0976*exp(-89.22*Delta) ...
                    - 0.0325*exp(-228*Delta) - 0.0509*exp(-125.9*Delta^(2./3.));
end

end