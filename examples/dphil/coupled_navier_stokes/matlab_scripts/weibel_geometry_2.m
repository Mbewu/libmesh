% plots using weibel geometry
% SI UNITS


%close all

% what plots to do
general_per_generation = true;
stk_per_generation = false;
deposition_per_generation = true;
v_d_ratio_per_generation = false;


rho = 1.2;
mu = 1.2*1.7e-5;

% want to be able to specify different geometries
% in cm and radians
%[diameter, branch_length, branching_angle, gravity_angle,generations] = weibel_data()
[diameter, branch_length, branching_angle, gravity_angle,generations] = james_symmetric_data()

% generation 0 should be the trachea
%generations= 0:23;  % generation 0 is the trachea 

% dependent stuff
num_branches = 2.^generations

area = pi*(diameter/2).^2;% pi*r^2 - cm^2
total_area = area.*num_branches; %total cross sectional area at each generation

inflow_rate = [10 20 40 120] * 1000/60; %cm^3/s 20lpm

branch_flow_rate_10 = inflow_rate(1)./num_branches; %flow_rate through each branch of a generation
branch_velocity_10 = branch_flow_rate_10./area;    %mean velocity in each branch
branch_reynolds_10 = rho*(branch_velocity_10*1e-2).*(diameter*1e-2)/mu;

branch_flow_rate_20 = inflow_rate(2)./num_branches; %flow_rate through each branch of a generation
branch_velocity_20 = branch_flow_rate_20./area;    %mean velocity in each branch
branch_reynolds_20 = rho*(branch_velocity_20*1e-2).*(diameter*1e-2)/mu;

branch_flow_rate_40 = inflow_rate(3)./num_branches; %flow_rate through each branch of a generation
branch_velocity_40 = branch_flow_rate_40./area;    %mean velocity in each branch
branch_reynolds_40 = rho*(branch_velocity_40*1e-2).*(diameter*1e-2)/mu;

branch_flow_rate_120 = inflow_rate(4)./num_branches; %flow_rate through each branch of a generation
branch_velocity_120 = branch_flow_rate_120./area;    %mean velocity in each branch
branch_reynolds_120 = rho*(branch_velocity_120*1e-2).*(diameter*1e-2)/mu;

if(general_per_generation)
    figure
    subplot(1,4,1)
    plot(generations,diameter);
    ylabel('Diameter (cm)');
    xlabel('Generation');
    xlim([0,23])
    legend
    subplot(1,4,2)
    hold on
    plot(generations,branch_flow_rate_10);
    plot(generations,branch_flow_rate_20);
    plot(generations,branch_flow_rate_40);
    plot(generations,branch_flow_rate_120);
    ylabel('flow_rate (cm^3/s)');
    xlabel('Generation');
    legend('10lpm - rest','20lpm - low activity','40lpm - light exertion','120lpm - heavy exertion');
    xlim([0,23])
    subplot(1,4,3)
    hold on
    plot(generations,branch_velocity_10);
    plot(generations,branch_velocity_20);
    plot(generations,branch_velocity_40);
    plot(generations,branch_velocity_120);
    legend('10lpm - rest','20lpm - low activity','40lpm - light exertion','120lpm - heavy exertion');
    ylabel('Mean Velocity (cm/s)');
    xlabel('Generation');
    xlim([0,23])
    subplot(1,4,4)
    hold on
    plot(generations,branch_reynolds_10);
    plot(generations,branch_reynolds_20);
    plot(generations,branch_reynolds_40);
    plot(generations,branch_reynolds_120);
    legend('10lpm - rest','20lpm - low activity','40lpm - light exertion','120lpm - heavy exertion');
    ylabel('Re');
    xlabel('Generation');
    xlim([0,23])
end


%%  stokes number SI units
rho_p = 1200;
lambda = 68.e-9
d_p = (5:11)*1e-6;
C = 1. + 2.*lambda./d_p.*(1.257 + 0.4*exp(-1.1*(d_p./(2*lambda))));

if(stk_per_generation)
    Stk_5 = rho_p*d_p(1)^2.*(branch_velocity_20(2:end)*1e-2)*C(1)./(18*mu*(diameter(2:end)*1e-2));
    Stk_6 = rho_p*d_p(2)^2.*(branch_velocity_20(2:end)*1e-2)*C(2)./(18*mu*(diameter(2:end)*1e-2));
    Stk_7 = rho_p*d_p(3)^2.*(branch_velocity_20(2:end)*1e-2)*C(3)./(18*mu*(diameter(2:end)*1e-2));
    Stk_8 = rho_p*d_p(4)^2.*(branch_velocity_20(2:end)*1e-2)*C(4)./(18*mu*(diameter(2:end)*1e-2));
    Stk_9 = rho_p*d_p(5)^2.*(branch_velocity_20(2:end)*1e-2)*C(5)./(18*mu*(diameter(2:end)*1e-2));
    Stk_10 = rho_p*d_p(6)^2.*(branch_velocity_20(2:end)*1e-2)*C(6)./(18*mu*(diameter(2:end)*1e-2));
    Stk_11 = rho_p*d_p(7)^2.*(branch_velocity_20(2:end)*1e-2)*C(7)./(18*mu*(diameter(2:end)*1e-2));

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
    title('Inflow flow rate 20lpm')
end



%% comparison of impaction and sedimentation

% stokes number SI units
rho_p = 1200;
lambda = 68.e-9
d_p = 1*1e-6;
C = 1. + 2.*lambda./d_p.*(1.257 + 0.4*exp(-1.1*(d_p./(2*lambda))));
diffusivity = 2.8e-11;  % for d = 1e-6m, tsuda 2013


stacked_bar_plot = true;
if(deposition_per_generation)
    
    St_10 = rho_p*d_p^2.*(branch_velocity_10*1e-2)*C./(18*mu*(diameter.*1e-2));

    p_imp_yeh_10 = yeh_impaction(branching_angle,St_10)
    p_sed_wang_10 = wang_sedimentation(St_10,gravity_angle,branch_length*1e-2,branch_velocity_10*1e-2)
    p_dif_ingham_10 = ingham_diffusion(branch_length*1e-2,branch_flow_rate_10*1e-6,diffusivity)
    
    St_20 = rho_p*d_p^2.*(branch_velocity_20*1e-2)*C./(18*mu*(diameter.*1e-2));

    p_imp_yeh_20 = yeh_impaction(branching_angle,St_20)
    p_sed_wang_20 = wang_sedimentation(St_20,gravity_angle,branch_length*1e-2,branch_velocity_20*1e-2)
    p_dif_ingham_20 = ingham_diffusion(branch_length*1e-2,branch_flow_rate_20*1e-6,diffusivity)

    St_40 = rho_p*d_p^2.*(branch_velocity_40*1e-2)*C./(18*mu*(diameter.*1e-2));

    p_imp_yeh_40 = yeh_impaction(branching_angle,St_40)
    p_sed_wang_40 = wang_sedimentation(St_40,gravity_angle,branch_length*1e-2,branch_velocity_40*1e-2)
    p_dif_ingham_40 = ingham_diffusion(branch_length*1e-2,branch_flow_rate_40*1e-6,diffusivity)

    St_120 = rho_p*d_p^2.*(branch_velocity_120*1e-2)*C./(18*mu*(diameter.*1e-2));

    p_imp_yeh_120 = yeh_impaction(branching_angle,St_120)
    p_sed_wang_120 = wang_sedimentation(St_120,gravity_angle,branch_length*1e-2,branch_velocity_120*1e-2)
    p_dif_ingham_120 = ingham_diffusion(branch_length*1e-2,branch_flow_rate_120*1e-6,diffusivity)

    figure
    subplot(1,4,1)
    plot(generations(1:15),p_imp_yeh_10(1:15))
    hold on
    plot(generations(1:15),p_sed_wang_10(1:15))
    hold on
    plot(generations(1:15),p_dif_ingham_10(1:15))
    title('flow 10lpm, particle size 2um')
    xlabel('generation')
    ylabel('deposition probability')
    legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
    xlim([0 23])
    subplot(1,4,2)
    plot(generations(1:15),p_imp_yeh_20(1:15))
    hold on
    plot(generations(1:15),p_sed_wang_20(1:15))
    hold on
    plot(generations(1:15),p_dif_ingham_20(1:15))
    title('flow 20lpm, particle size 2um')
    xlabel('generation')
    ylabel('deposition probability')
    legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
    xlim([0 23])
    subplot(1,4,3)
    plot(generations(1:15),p_imp_yeh_40(1:15))
    hold on
    plot(generations(1:15),p_sed_wang_40(1:15))
    hold on
    plot(generations(1:15),p_dif_ingham_40(1:15))
    title('flow 40lpm, particle size 2um')
    xlabel('generation')
    ylabel('deposition probability')
    legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
    xlim([0 23])
    subplot(1,4,4)
    plot(generations(1:15),p_imp_yeh_120(1:15))
    hold on
    plot(generations(1:15),p_sed_wang_120(1:15))
    hold on
    plot(generations(1:15),p_dif_ingham_120(1:15))
    title('flow 120lpm, particle size 1um')
    xlabel('generation')
    ylabel('deposition probability')
    legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
    xlim([0 23])
    
    if(stacked_bar_plot)
        
        p_matrix_10 = [p_imp_yeh_10;p_sed_wang_10;p_dif_ingham_10]'
        [rows, columns] = size(p_matrix_10);
        denom = repmat(sum(p_matrix_10,2),[1,columns])
        p_matrix_10_norm = p_matrix_10./denom;
        
        p_matrix_20 = [p_imp_yeh_20;p_sed_wang_20;p_dif_ingham_20]'
        [rows, columns] = size(p_matrix_20);
        denom = repmat(sum(p_matrix_20,2),[1,columns])
        p_matrix_20_norm = p_matrix_20./denom;
        
        
        p_matrix_40 = [p_imp_yeh_40;p_sed_wang_40;p_dif_ingham_40]'
        [rows, columns] = size(p_matrix_40);
        denom = repmat(sum(p_matrix_40,2),[1,columns])
        p_matrix_40_norm = p_matrix_40./denom;
        
        p_matrix_120 = [p_imp_yeh_120;p_sed_wang_120;p_dif_ingham_120]'
        [rows, columns] = size(p_matrix_120);
        denom = repmat(sum(p_matrix_120,2),[1,columns])
        p_matrix_120_norm = p_matrix_120./denom;
        
        figure
        subplot(1,4,1)
        bar(generations,p_matrix_10_norm,'stacked')
        title('flow 10lpm, particle size 1um')
        xlabel('generation')
        ylabel('deposition proportion')
        legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
        xlim([0 23])
        ylim([0 1.2])
        subplot(1,4,2)
        bar(generations,p_matrix_20_norm,'stacked')
        title('flow 20lpm, particle size 1um')
        xlabel('generation')
        ylabel('deposition proportion')
        legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
        xlim([0 23])
        ylim([0 1.2])
        subplot(1,4,3)
        bar(generations,p_matrix_40_norm,'stacked')
        title('flow 40lpm, particle size 1um')
        xlabel('generation')
        ylabel('deposition proportion')
        legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
        xlim([0 23])
        ylim([0 1.2])
        subplot(1,4,4)
        bar(generations,p_matrix_120_norm,'stacked')
        title('flow 120lpm, particle size 1um')
        xlabel('generation')
        ylabel('deposition proportion')
        legend('impaction (yeh)','sedimentation (wang)','diffusion (ingham)')
        xlim([0 23])
        ylim([0 1.2])
    end
end



%% velocity magnitude diameter ratio

if(v_d_ratio_per_generation)
    plot(generations,branch_velocity_20./diameter,'*');
    ylabel('u_0 / D (s^{{-}1})');
    xlabel('Generation');
    xlim([0,23])
    title('Inflow flow rate 20lpm')
    set(gca,'fontsize', 18);
end
