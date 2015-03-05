%plot some pressure stuff

flow = 60*[0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.75 1.0 1.25 1.5 1.75 2.0];
flow_string = ['0.01'; '0.05'; '0.1 '; '0.2 '; '0.3 '; '0.4 '; '0.5 '; '0.75'; '1.0 '; '1.25'; '1.5 '; '1.75'; '2.0 '];
celldata = cellstr(flow_string);
poiseuille = flow;
pedley = flow;
ertbruggen = flow;
reynolds = flow;

for i=1:length(flow)
    A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/pressure_flow_study/poiseuille/',celldata{i},'/out.dat'));
    poiseuille(i) = A.data(10,5)*0.01019;
end

for i=1:length(flow)
    A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/pressure_flow_study/pedley/',celldata{i},'/out.dat'));
    pedley(i) = A.data(10,5)*0.01019;
end

for i=1:length(flow)
    A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/pressure_flow_study/ertbruggen/',celldata{i},'/out.dat'));
    ertbruggen(i) = A.data(10,5)*0.01019;
end


for i=1:length(flow)
    A = importdata(strcat('~/libmesh-0.9.3/examples/dphil/coupled_navier_stokes/results/pressure_flow_study/reynolds/',celldata{i},'/out.dat'));
    reynolds(i) = A.data(10,5)*0.01019;
end

poiseuille
pedley
ertbruggen
reynolds

figure
hold on
plot(flow,poiseuille,'r*-');
plot(flow,pedley,'b*-');
plot(flow,ertbruggen,'m*-');
plot(flow,reynolds,'g*-');
legend('poiseuille','pedley','ertbruggen','reynolds');
xlabel('flow rate (L/min)')
ylabel('pressure drop from trachea (cmH2O)')
title('pressure-flow on tree APLE0081560')
