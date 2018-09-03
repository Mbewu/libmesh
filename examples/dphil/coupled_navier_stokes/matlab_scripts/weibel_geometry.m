
generations= 0:23;  % generation 0 is the trachea 
rho = 1.2;
mu = 1.2*1.7e-5;

num_branches = 2.^generations
diameter = [1.8 1.22 0.83 0.56 0.45 0.35 0.28 0.23 0.186 0.154 0.13 0.109 ...
    0.095 0.082 0.074 0.066 0.06 0.054 0.05 0.047 0.045 0.043 0.041 0.041]; %cm
branch_length = [12 4.76 1.9 0.76 1.27 1.07 0.9 0.76 0.64 0.54 0.46 0.39 0.33 0.27 ...
    0.23 0.2 0.165 0.141 0.117 0.099 0.083 0.07 0.059 0.05]; % cm
area = pi*(diameter/2).^2;% pi*r^2 - cm^2
total_area = area.*num_branches; %total cross sectional area at each generation

%yeh 1980 angles
branching_angle = pi/180*[0 33 34 22 20 18 19 22 28 22 33 34 37 39 39 51 45 45 45 45 45 45 45 45];
gravity_angle =   pi/180*[0 20 31 43 39 39 40 36 39 45 43 45 45 60 60 60 60 60 60 60 60 60 60 60];

influx = [20 40 120] * 1000/60; %cm^3/s 20lpm
branch_flux_20 = influx(1)./num_branches; %flux through each branch of a generation
branch_velocity_20 = branch_flux_20./area;    %mean velocity in each branch
branch_reynolds_20 = rho*(branch_velocity_20*1e-2).*(diameter*1e-2)/mu;

branch_flux_40 = influx(2)./num_branches; %flux through each branch of a generation
branch_velocity_40 = branch_flux_40./area;    %mean velocity in each branch
branch_reynolds_40 = rho*(branch_velocity_40*1e-2).*(diameter*1e-2)/mu;

branch_flux_120 = influx(3)./num_branches; %flux through each branch of a generation
branch_velocity_120 = branch_flux_120./area;    %mean velocity in each branch
branch_reynolds_120 = rho*(branch_velocity_120*1e-2).*(diameter*1e-2)/mu;

close all
subplot(1,4,1)
plot(generations,diameter);
ylabel('Diameter (cm)');
xlabel('Generation');
xlim([0,23])
legend
subplot(1,4,2)
hold on
plot(generations,branch_flux_20);
plot(generations,branch_flux_40);
plot(generations,branch_flux_120);
ylabel('Flux (cm^3/s)');
xlabel('Generation');
legend('20l/m - low activity','40l/m - light exertion','120l/m - heavy exertion');
xlim([0,23])
subplot(1,4,3)
hold on
plot(generations,branch_velocity_20);
plot(generations,branch_velocity_40);
plot(generations,branch_velocity_120);
legend('20l/m - low activity','40l/m - light exertion','120l/m - heavy exertion');
ylabel('Mean Velocity (cm/s)');
xlabel('Generation');
xlim([0,23])
subplot(1,4,4)
hold on
plot(generations,branch_reynolds_20);
plot(generations,branch_reynolds_40);
plot(generations,branch_reynolds_120);
legend('20l/m - low activity','40l/m - light exertion','120l/m - heavy exertion');
ylabel('Re');
xlabel('Generation');
xlim([0,23])

% stokes number SI units
rho_p = 1200;
lambda = 68.e-9
d_p = (5:11)*1e-6;
C = 1. + 2.*lambda./d_p.*(1.257 + 0.4*exp(-1.1*(d_p./(2*lambda))));

Stk_5 = rho_p*d_p(1)^2.*(branch_velocity(2:end)*1e-2)*C(1)./(18*mu*(diameter(2:end)*1e-2));
Stk_6 = rho_p*d_p(2)^2.*(branch_velocity(2:end)*1e-2)*C(2)./(18*mu*(diameter(2:end)*1e-2));
Stk_7 = rho_p*d_p(3)^2.*(branch_velocity(2:end)*1e-2)*C(3)./(18*mu*(diameter(2:end)*1e-2));
Stk_8 = rho_p*d_p(4)^2.*(branch_velocity(2:end)*1e-2)*C(4)./(18*mu*(diameter(2:end)*1e-2));
Stk_9 = rho_p*d_p(5)^2.*(branch_velocity(2:end)*1e-2)*C(5)./(18*mu*(diameter(2:end)*1e-2));
Stk_10 = rho_p*d_p(6)^2.*(branch_velocity(2:end)*1e-2)*C(6)./(18*mu*(diameter(2:end)*1e-2));
Stk_11 = rho_p*d_p(7)^2.*(branch_velocity(2:end)*1e-2)*C(7)./(18*mu*(diameter(2:end)*1e-2));

figure
plot(generations(2:end),Stk_5)
hold on
plot(generations(2:end),Stk_6)
plot(generations(2:end),Stk_7)
plot(generations(2:end),Stk_8)
plot(generations(2:end),Stk_9)
plot(generations(2:end),Stk_10)
plot(generations(2:end),Stk_11)
legend('5um','6um','7um','8um','9um','10um','11um')
ylabel('Stokes number')



%comparison of impaction and sedimentation

% stokes number SI units
rho_p = 1200;
lambda = 68.e-9
d_p = 10*1e-6;
C = 1. + 2.*lambda./d_p.*(1.257 + 0.4*exp(-1.1*(d_p./(2*lambda))));

St_20 = rho_p*d_p^2.*(branch_velocity_20*1e-2)*C./(18*mu*(diameter.*1e-2));

p_imp_yeh_20 = yeh_impaction(branching_angle,St_20)
p_sed_wang_20 = wang_sedimentation(St_20,gravity_angle,branch_length*1e-2,branch_velocity_20*1e-2)

St_40 = rho_p*d_p^2.*(branch_velocity_40*1e-2)*C./(18*mu*(diameter.*1e-2));

p_imp_yeh_40 = yeh_impaction(branching_angle,St_40)
p_sed_wang_40 = wang_sedimentation(St_40,gravity_angle,branch_length*1e-2,branch_velocity_40*1e-2)

St_120 = rho_p*d_p^2.*(branch_velocity_120*1e-2)*C./(18*mu*(diameter.*1e-2));

p_imp_yeh_120 = yeh_impaction(branching_angle,St_120)
p_sed_wang_120 = wang_sedimentation(St_120,gravity_angle,branch_length*1e-2,branch_velocity_120*1e-2)



figure
subplot(1,3,1)
plot(generations,p_imp_yeh_20)
hold on
plot(generations,p_sed_wang_20)
title('comparison of deposition mechanisms - flow 20lpm')
legend('impaction (yeh)','sedimentation (wang)')
xlim([0 23])
subplot(1,3,2)
plot(generations,p_imp_yeh_40)
hold on
plot(generations,p_sed_wang_40)
title('comparison of deposition mechanisms - flow 40lpm')
legend('impaction (yeh)','sedimentation (wang)')
xlim([0 23])
subplot(1,3,3)
plot(generations,p_imp_yeh_120)
hold on
plot(generations,p_sed_wang_120)
title('comparison of deposition mechanisms - flow 120lpm')
legend('impaction (yeh)','sedimentation (wang)')
xlim([0 23])
