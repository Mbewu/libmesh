% most taken from lambert

pi = 3.142

g = 9.81;   %m/s^2
C = 1.0;
rho_p = 1200;
r_p = 30.e-6; %2.5,10,30 um
L = 2.0e-2; %1cm
phi = 0.;
mu = 1.2*1.7e-5;    %m^2/s
R = 0.5e-2; % 0.5cm
v = 24.80e-2;   %cm/s
pi = 3.142;

r_p = linspace(0e-6,20e-6,11)


probability_wang = 1 - exp(-4*g*C*rho_p.*r_p.*r_p*L*cos(phi)/(9*pi*mu*R*v))

close all
plot(r_p,probability_wang)

%dt 0.001
prob_exp_1 = [0 0.026 0.096 0.148 0.214 0.265 0.316 0.366 0.401 0.461 0.507]



 
hold on
plot(r_p,prob_exp_1)
%plot(r_p,prob_exp_1)
%plot(r_p,prob_exp_2)
%plot(r_p,prob_exp_3)
%plot(r_p,prob_exp_4)
xlabel('particle radius');
ylabel('deposition probability');
legend('Wang 0D','3D mean 0.001')
ax = gca;
ax.XAxis.Exponent = -6;
