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
r_p = logspace(log10(0.1e-6),log10(50e-6),1001);

St_100 = C*rho_p.*r_p.*r_p*v_100/(9*mu*R);
St_200 = C*rho_p.*r_p.*r_p*v_200/(9*mu*R);
St_400 = C*rho_p.*r_p.*r_p*v_400/(9*mu*R);
St_1000 = C*rho_p.*r_p.*r_p*v_1000/(9*mu*R);

p_imp_zhang_100 = zhang_variable_radius(r_p,v_100);
p_imp_zhang_200 = zhang_variable_radius(r_p,v_200);
p_imp_zhang_400 = zhang_variable_radius(r_p,v_400);
p_imp_zhang_1000 = zhang_variable_radius(r_p,v_1000);



%%


% what particle radius corrsponds to stk=1
r_p_stk_1 = sqrt(9*mu*R/(C*rho_p)); 
% what particle radius corrsponds to stk=0.001
r_p_stk_0001 = sqrt(0.001*9*mu*R/(C*rho_p)); 



close all
figure
hold on
plot(r_p,p_imp_zhang_100)
plot(r_p,p_imp_zhang_200)
plot(r_p,p_imp_zhang_400)
plot(r_p,p_imp_zhang_1000)
xlabel('particle radius');
ylabel('deposition probability');
title('impaction models curved 45 angle')
legend('zhang, vmax=100','zhang, vmax=200','zhang, vmax=400','zhang, vmax=1000')
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
plot(St_100,p_imp_zhang_100)
hold on
plot(St_200,p_imp_zhang_200)
plot(St_400,p_imp_zhang_400)
plot(St_1000,p_imp_zhang_1000)
xlabel('Stk');
ylabel('deposition probability');
title('impaction models curved 45 angle')
legend('zhang, vmax=100','zhang, vmax=200','zhang, vmax=400','zhang, vmax=1000')
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
plot(St_100,p_imp_zhang_100)
hold on
plot(St_200,p_imp_zhang_200)
plot(St_400,p_imp_zhang_400)
plot(St_1000,p_imp_zhang_1000)
xlabel('Stk');
ylabel('deposition probability');
title('impaction models curved 45 angle')
legend('zhang, vmax=100','zhang, vmax=200','zhang, vmax=400','zhang, vmax=1000')
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