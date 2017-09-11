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
v = 50e-2;%99.20e-2 * 4;   %cm/s average velocity
pi = 3.142;
theta = pi/4;

%r_p = linspace(0e-6,20e-6,11)
%r_p = linspace(0e-6,50e-6,1001)
r_p = logspace(log10(0.1e-6),log10(50e-6),1001);

St = C*rho_p.*r_p.*r_p*v/(9*mu*R);
Re = rho*v*2*R/mu;


%yeh impaction
prob_imp_yeh = 0;

theta*St
%prob_imp_yeh(1) = 0.;
%for i=2:length(r_p)
for i=1:length(r_p)
    if theta*St(i) < 1
        prob_imp_yeh(i) = 1 - 2/pi.*acos(theta*St(i)) + 1./pi*sin(2.*acos(theta*St(i)));
    else
        prob_imp_yeh(i) = 1;
    end
end

prob_imp_yeh

%zhang impaction - parabolic
for i=1:length(St)
    if(St(i) < 0.04)
        prob_imp_zhang(i) = 0.000654*exp(55.7*St(i)^0.954)*Re^(1/3)*sin(theta);
    else
        prob_imp_zhang(i) = (0.19 - 0.193*exp(-9.5*St(i)^1.565))*Re^(1/3)*sin(theta);
    end
end

%cai impaction
prob_imp_cai = 0;
R_d = R;
f_0 = pi*(1-1/4*(R_d/R)^2)-(4/3*(15/16*pi - 2)*(R_d/R)^2)*(cos(theta))^2;
f_1 = 1 - 1/3*(R_d/R)^2 + (pi - 11/3)*(R_d/R)^2*(cos(theta))^2 ...
    -1/3*(R_d/R)^2*sin(theta) + (2/3 - pi/8)*(R_d/R)^4*(cos(theta))^2 ...
    +1/5*(R_d/R)^4*sin(theta)^2 + (6 - 15/8*pi)*(R_d/R)^4*(cos(theta))^4 ...
    +(7/15 - pi/8)*(R_d/R)^4*(sin(theta))^2*(cos(theta))^2;
G = 8*sin(theta)*f_1/((R_d/R)*f_0);
prob_imp_cai = G*St;


%experimental kim
St_kim = [0.05 0.14 0.27];
prob_imp_kim = [0.063 0.375 0.684]

%computational robinson 1,2
prob_imp_rob_1 = [0.11 0.30 0.54]
prob_imp_rob_2 = [0.04 0.21 0.49]

close all

%% SIMULATION RESULTS

% particle radius
%r_p_sim = [0 5 10 15 20 25 30 35 40 45 50]*1e-6;
r_p_sim = [0 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100]*1e-6;

St_sim = C*rho_p.*r_p_sim.*r_p_sim*v/(9*mu*R);

%dt 0.001
prob_exp_1 = r_p_sim*0;
prob_exp_no_bdy = r_p_sim*0;


% 0,10,20,...,100
% curved - no bdy - 45 - 200cm/s
prob_exp_no_bdy(1:11) = [0 0.006 0.043 0.091 0.339 0.542 0.627 0.709 0.737 0.793 0.838]

% curved - bdy layer - 45 - 200cm/s
prob_exp_1(1:11) = [0 0.011 0.052 0.137 0.399 0.526 0.650 0.729 0.798 0.801 0.851]


% curved - no bdy - 45 - 100cm/s
%prob_exp_1(1:11) = [0 0.005 0.036 0.099 0.189 0.428 0.573 0.644 0.693 0.703 0.787]

% curved - bdy layer - 45 - 100cm/s
%prob_exp_1(1:11) = [0 0.054 0.112 0.154 0.323 0.458 0.577 0.661 0.735 0.807 0.819]

% straight - bdy - 45 - 100cm/s
% 0 1 2 3 4 5 6 7 8 9 10 20 25 30 40 50 60 70 80 90
%prob_exp_1 = [0 0.014 0.014 0.009 0.019 0.016 0.020 0.025 0.017 0.018 0.026 0.081 0.258 0.432 0.536 0.648 0.731 0.777 0.817 0.836]
% 0 1 10 20 25 30 40 50 60 70 80 90 100
prob_exp_100 = [0 0.014 0.026 0.081 0.258 0.450 0.573 0.691 0.768 0.819 0.864 0.884];
prob_exp_200 = [0 0.014 0.026 0.081 0.258 0.450 0.573 0.691 0.768 0.819 0.864 0.884];
prob_exp_400 = [0 0.014 0.026 0.081 0.258 0.450 0.573 0.691 0.768 0.819 0.864 0.884];
prob_exp_1000 = [0 0.014 0.026 0.081 0.258 0.450 0.573 0.691 0.768 0.819 0.864 0.884];



%%


% what particle radius corrsponds to stk=1
r_p_stk_1 = sqrt(9*mu*R/(C*rho_p)); 
% what particle radius corrsponds to stk=0.001
r_p_stk_0001 = sqrt(0.001*9*mu*R/(C*rho_p)); 



 
hold on
plot(r_p,prob_imp_zhang)
plot(r_p,prob_imp_yeh)
plot(r_p,prob_imp_cai)
plot(r_p_sim,prob_exp_1,'*-')
plot(r_p_sim,prob_exp_no_bdy,'*-')
%plot(r_p,prob_exp_1)
%plot(r_p,prob_exp_2)
%plot(r_p,prob_exp_3)
%plot(r_p,prob_exp_4)
xlabel('particle radius');
ylabel('deposition probability');
title('impaction models curved 45 angle')
legend('model - zhang','model - yeh','model - cai','mysim','mysim no bdy')
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
plot(St,prob_imp_zhang)
hold on
plot(St,prob_imp_yeh)
plot(St,prob_imp_cai)
plot(St_kim,prob_imp_kim,'*')
plot(St_kim,prob_imp_rob_1,'+')
plot(St_kim,prob_imp_rob_2,'o')
plot(St_sim,prob_exp_1,'*-')
plot(St_sim,prob_exp_no_bdy,'*-')
xlabel('Stk');
ylabel('deposition probability');
title('impaction models curved 45 angle')
legend('model - zhang','model - yeh','model - cai','exp - kim','sim - rob 1','sim - rob 2','mysim','mysim no bdy')
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
plot(St,prob_imp_zhang)
hold on
plot(St,prob_imp_yeh)
plot(St,prob_imp_cai)
plot(St_kim,prob_imp_kim,'*')
plot(St_kim,prob_imp_rob_1,'+')
plot(St_kim,prob_imp_rob_2,'o')
plot(St_sim,prob_exp_1,'*-')
plot(St_sim,prob_exp_no_bdy,'*-')
xlabel('Stk');
ylabel('deposition probability');
title('impaction models curved 45 angle')
legend('model - zhang','model - yeh','model - cai','exp - kim','sim - rob 1','sim - rob 2','mysim','mysim no bdy')
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