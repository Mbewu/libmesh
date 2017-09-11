% most taken from lambert

log_plot = false;

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
v_100 = 50e-2;%99.20e-2 * 4;   %cm/s average velocity
v_200 = 100e-2;%99.20e-2 * 4;   %cm/s average velocity
v_400 = 200e-2;%99.20e-2 * 4;   %cm/s average velocity
v_1000 = 500e-2;%99.20e-2 * 4;   %cm/s average velocity
pi = 3.142;
theta = pi/4;

%r_p = linspace(0e-6,20e-6,11)
%r_p = linspace(0e-6,50e-6,1001)
r_p = [0 1 10 20 30 40 50 60 70 80 90 100]*0.5e-6;

St_100 = C*rho_p.*r_p.*r_p*v_100/(9*mu*R);
St_200 = C*rho_p.*r_p.*r_p*v_200/(9*mu*R);
St_400 = C*rho_p.*r_p.*r_p*v_400/(9*mu*R);
St_1000 = C*rho_p.*r_p.*r_p*v_1000/(9*mu*R);

% straight - bdy - 45 - 100cm/s
% 0 1 2 3 4 5 6 7 8 9 10 20 25 30 40 50 60 70 80 90
%prob_exp_1 = [0 0.014 0.014 0.009 0.019 0.016 0.020 0.025 0.017 0.018 0.026 0.081 0.258 0.432 0.536 0.648 0.731 0.777 0.817 0.836]
% 0 1 10 20 25 30 40 50 60 70 80 90 100
prob_exp_100 = [0 0.014 0.026 0.081 0.258 0.450 0.573 0.691 0.768 0.819 0.864 0.884];
prob_exp_200 = [0 0.018 0.023 0.100 0.330 0.501 0.641 0.737 0.796 0.852 0.892 0.917];
prob_exp_400 = [0 0.010 0.028 0.130 0.336 0.533 0.662 0.763 0.839 0.851 0.903  0.904];
prob_exp_1000 = [0 0.015 0.038 0.173 0.419 0.608 0.733 0.763 0.824 0.831 0.852 0.861];



%%


% what particle radius corrsponds to stk=1
r_p_stk_1 = sqrt(9*mu*R/(C*rho_p)); 
% what particle radius corrsponds to stk=0.001
r_p_stk_0001 = sqrt(0.001*9*mu*R/(C*rho_p)); 



close all
figure
hold on
plot(r_p,prob_exp_100,'*-')
plot(r_p,prob_exp_200,'*-')
plot(r_p,prob_exp_400,'*-')
plot(r_p,prob_exp_1000,'*-')
xlabel('particle radius');
ylabel('deposition probability');
title('impaction mysim curved 45 angle')
legend('mysim, vmax=100','mysim, vmax=200','mysim, vmax=400','mysim, vmax=1000')
ax = gca;
ax.XAxis.Exponent = -6;
ylim([0 1.5])
% make line at stk = 1
yL = get(gca,'YLim');
line([r_p_stk_1 r_p_stk_1],yL,'Color','k');

% log plot
if(log_plot)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([r_p_stk_0001 max(r_p)])
end


%impaction
figure
plot(St_100,prob_exp_100,'*-')
hold on
plot(St_200,prob_exp_200,'*-')
plot(St_400,prob_exp_400,'*-')
plot(St_1000,prob_exp_1000,'*-')
xlabel('Stk');
ylabel('deposition probability');
title('impaction mysim curved 45 angle')
legend('mysim, vmax=100','mysim, vmax=200','mysim, vmax=400','mysim, vmax=1000')
ylim([0 1.5])
% make line at stk = 1
yL = get(gca,'YLim');
line([1 1],yL,'Color','k');

% log plot
if(log_plot)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1e-3 max(St)])
end


%impaction
figure
plot(St_100,prob_exp_100,'*-')
hold on
plot(St_200,prob_exp_200,'*-')
plot(St_400,prob_exp_400,'*-')
plot(St_1000,prob_exp_1000,'*-')
xlabel('Stk');
ylabel('deposition probability');
title('impaction mysim curved 45 angle')
legend('mysim, vmax=100','mysim, vmax=200','mysim, vmax=400','mysim, vmax=1000')
xlim([0 0.3])
ylim([0 1.0])
xL = get(gca,'XLim');
line([1 1],xL,'Color','k');

% log plot
if(log_plot)
    set(gca,'xscale','log');
    set(gca,'yscale','log');
    xlim([1e-3 0.3])
end