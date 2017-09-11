% most taken from lambert

pi = 3.142

g = 9.81;   %m/s^2
C = 1.0;
rho = 1.2;
rho_p = 1200;
r_p = 30.e-6; %2.5,10,30 um
L = 5.3e-2; %1cm
phi = 0.;
mu = 1.2*1.7e-5;    %m^2/s
R = 0.5e-2; % 0.5cm
R_d = 0.5e-2;
v = 24.80e-2;   %cm/s
pi = 3.142;
theta = pi/4;

%r_p = linspace(0e-6,20e-6,11)
r_p = linspace(0e-6,50e-6,11)

St = C*rho_p.*r_p.*r_p*v/(9*mu*R)
Re = rho*v*2*R/mu


%yeh impaction
prob_imp_yeh = 0;
if theta*St < 1
    prob_imp_yeh = 1 - 2/pi.*acos(theta*St) + 1./pi*sin(2.*acos(theta*St))
else
    prob_imp_yeh = 1
end

%zhang impaction - parabolic
prob_imp_zhang = 0;
if(St < 0.04)
    prob_imp_zhang = 0.000654*exp(55.7*St.^0.954)*Re^(1/3)*sin(theta)
else
    prob_imp_zhang = (0.19 - 0.193*exp(-9.5*St.^1.565))*Re^(1/3)*sin(theta)
end

%cai impaction
prob_imp_cai = 0;
f_0 = pi*(1-1/4*(R_d/R)^2)-(4/3*(15/16*pi - 2)*(R_d/R)^2)*(cos(theta))^2;
f_1 = 1 - 1/3*(R_d/R)^2 + (pi - 11/3)*(R_d/R)^2*(cos(theta))^2 ...
    -1/3*(R_d/R)^2*sin(theta) + (2/3 - pi/8)*(R_d/R)^4*(cos(theta))^2 ...
    +1/5*(R_d/R)^4*sin(theta)^2 + (6 - 15/8*pi)*(R_d/R)^4*(cos(theta))^4 ...
    +(7/15 - pi/8)*(R_d/R)^4*(sin(theta))^2*(cos(theta))^2;
G = 8*sin(theta)*f_1/((R_d/R)*f_0);
prob_imp_cai = G*St


%experimental kim
St_kim = [0.05 0.14 0.27];
prob_imp_kim = [0.063 0.375 0.684]

%computational robinson 1,2
prob_imp_rob_1 = [0.11 0.30 0.54]
prob_imp_rob_2 = [0.04 0.21 0.49]

%dt 0.001
prob_exp_1 = r_p*0;
prob_exp_1(1:11) = [0 0.0 0.009 0.034 0.059 0.078 0.150 0.294 0.363 0.457  0.540]




%close all
figure
hold on
plot(r_p,prob_imp_zhang)
plot(r_p,prob_imp_yeh)
plot(r_p,prob_imp_cai)
plot(r_p,prob_exp_1)
%plot(r_p,prob_exp_1)
%plot(r_p,prob_exp_2)
%plot(r_p,prob_exp_3)
%plot(r_p,prob_exp_4)
xlabel('particle radius');
ylabel('deposition probability');
legend('imp_zhang','imp_yeh','imp_cai','3D mean 0.001')
ax = gca;
ax.XAxis.Exponent = -6;
ylim([0 1.0])


%impaction
figure
plot(St,prob_imp_zhang)
hold on
plot(St,prob_imp_yeh)
plot(St,prob_imp_cai)
plot(St_kim,prob_imp_kim,'*')
plot(St_kim,prob_imp_rob_1,'+')
plot(St_kim,prob_imp_rob_2,'o')
plot(St,prob_exp_1)
xlabel('Stk');
ylabel('deposition probability');
legend('imp zhang','imp yeh','imp cai','imp kim','imp rob 1','imp rob 2','3D mean 0.001')
title('impaction models')
ylim([0 1.0])

%impaction
figure
plot(St,prob_imp_zhang)
hold on
plot(St,prob_imp_yeh)
plot(St,prob_imp_cai)
plot(St_kim,prob_imp_kim,'*')
plot(St_kim,prob_imp_rob_1,'+')
plot(St_kim,prob_imp_rob_2,'o')
plot(St,prob_exp_1)
xlabel('Stk');
ylabel('deposition probability');
legend('imp zhang','imp yeh','imp cai','imp kim','imp rob 1','imp rob 2','3D mean 0.001')
title('impaction models')
xlim([0 0.3])
ylim([0 1.0])