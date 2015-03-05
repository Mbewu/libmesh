%plot some pressure stuff

time_01 = 0:0.1:0.5;
time_001 = 0:0.01:0.5;
time_0002 = 0:0.002:0.5;
flow_string = ['0.01'; '0.05'; '0.1 '; '0.2 '; '0.3 '; '0.4 '; '0.5 '; '0.75'; '1.0 '; '1.25'; '1.5 '; '1.75'; '2.0 '];
celldata = cellstr(flow_string);

A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/coupling_study_new/dt_0.1/monolithic/out.dat'));
mono_time_01 = A.data(:,2);
monolithic_fluxes_3d = A.data(:,3:8);
monolithic_pressures_3d = A.data(:,9:14);
monolithic_fluxes_1d = A.data(:,15:19);
monolithic_pressures_1d = A.data(:,20:24);  

% 
% for i=1:length(flow)
%     A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/pressure_flow_study/ertbruggen/',celldata{i},'/out.dat'));
%     ertbruggen(i) = A.data(10,5)*0.01019;
% end
% 
% 
% for i=1:length(flow)
%     A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/pressure_flow_study/reynolds/',celldata{i},'/out.dat'));
%     reynolds(i) = A.data(10,5)*0.01019;
% end
% 
% poiseuille
% pedley
% ertbruggen
% reynolds

figure
plot(mono_time_01,monolithic_fluxes_3d(:,1));
title('fluxes');
xlabel('time');
ylabel('fluxes');
