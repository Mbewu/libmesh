% most taken from lambert

pi = 3.142

g = 9.81;   %m/s^2
C = 1.0;
rho_p = 1200;
r_p = 30.e-6; %2.5,10,30 um
L = 1.0e-2; %1cm
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
prob_exp_1 = [0 0.022 0.050 0.110 0.140 0.192 0.203 0.26 0.2698 0.3143 0.342685]
prob_exp_2 = [0 0.018 0.051 0.136 0.131 0.194 0.213 0.255 0.2605 0.3133 0.329]
prob_exp_3 = [0 0.014 0.058 0.098 0.131 0.197 0.206 0.248 0.2796 0.3 0.311]
prob_exp_4 = [0 0.013 0.055 0.092 0.114 0.189 0.191 0.2565 0.295 0.3049 0.323]

prob_exp = (prob_exp_1 + prob_exp_2 + prob_exp_3 + prob_exp_4)/4;

%dt 0.01
prob_exp_dt_01 = [0 0.152 0.181 0.239 0.264 0.298 0.334 0.340 0.354 0.382 0.394]

%dt 0.01
prob_exp_dt_0001 = [0 0.014 0.050 0.096 0.143 0.181 0.214 0.233 0.292 0.308 0.316]


 
hold on
plot(r_p,prob_exp)
plot(r_p,prob_exp_dt_01)
plot(r_p,prob_exp_dt_0001)
%plot(r_p,prob_exp_1)
%plot(r_p,prob_exp_2)
%plot(r_p,prob_exp_3)
%plot(r_p,prob_exp_4)
xlabel('particle radius');
ylabel('deposition probability');
legend('Wang 0D','3D mean 0.001','3D 0.01','3D 0.0001')
ax = gca;
ax.XAxis.Exponent = -6;
