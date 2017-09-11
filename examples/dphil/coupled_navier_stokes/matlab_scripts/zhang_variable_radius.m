% most taken from lambert
% returns the zhang deposition probability as a function of radius for a
% given set of other parameters

function [p_imp_zhang] = zhang_variable_radius(r_p,v)

pi = 3.142

g = 9.81;   %m/s^2
C = 1.0;
rho = 1.2;
rho_p = 1200;
L = 5.3e-2; %1cm
phi = 0.;
mu = 1.2*1.7e-5;    %m^2/s
R = 0.5e-2; % 0.5cm
pi = 3.142;
theta = pi/4;

St = C*rho_p.*r_p.*r_p*v/(9*mu*R);
Re = rho*v*2*R/mu;


%zhang impaction - parabolic
for i=1:length(St)
    if(St(i) < 0.04)
        p_imp_zhang(i) = 0.000654*exp(55.7*St(i)^0.954)*Re^(1/3)*sin(theta);
    else
        p_imp_zhang(i) = (0.19 - 0.193*exp(-9.5*St(i)^1.565))*Re^(1/3)*sin(theta);
    end
end
