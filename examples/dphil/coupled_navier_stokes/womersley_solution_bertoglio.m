clear
% womersley solution
axi = true;
dx = 0.0025;
dt = 0.002;
Re = 1./0.035; %10;
mu = 0.035; %3.5;
P = 4000;
R = 0.5; %0.5;
r = R:-dx:0;
real_r = [-r fliplr(r(1:end-1))];
if(axi)
    real_r = fliplr(r(1:end));
end
real_r_norm = real_r/R;

if(axi)
    r = fliplr(r(1:end));    
else
    r = [r fliplr(r(1:end-1))];    
end
u0 = -200.0/pi * (1 - besselj(0,(R*sqrt(4*pi/mu) * i^(1.5) * r/R)) ...
                    /besselj(0,(R*sqrt(4*pi/mu) * i^(1.5))));  %n = 1
%u0 = -200.0/pi * (1 - besselj(0,(R*sqrt(4*pi/mu) * i^0.5 * i^(1.5) * r/R)) ...
%                    /besselj(0,(R*sqrt(4*pi/mu) * i^0.5 * i^(1.5))));       %n = -1
t = 0:dt:1.0;
u = zeros(length(t),length(u0));
u = conj(cos(4*pi*t) + i*sin(4*pi*t))' * u0;   %n = 1
%u = (cos(4*pi*t) - i*sin(4*pi*t))' * u0;   n = -1
u = real(u);

max_mag = max(abs(max(max(real(u)))),abs(min(min(real(u)))));

%%%%%%%%%%%%%%%%%% plot solution
%u(1,:)
max(abs(u(1,:)))
max(abs(u(2,:)))
max(abs(u(3,:)))

figure
hold off
for j=1:length(t)
   plot([-fliplr(real_r) real_r(2:end)],[fliplr(real(u(j,:))) real(u(j,2:end))]);
   axis([-real_r(end) real_r(end) -max_mag max_mag])
   %plot(real(u0(:)));
   title(['t = ' num2str(t(j))])
   pause(0.1)
end

%%%%%%%%%%%%%%%%%% write file
% we want to write it with normalised R
% filename = 'womersley_solution_axi_coarse.dat';
% fid = fopen(filename,'w');
% fprintf(fid,'-1');
% for j=1:length(real_r_norm)
%     fprintf(fid,' %.7E',real_r_norm(j));
% end
% fprintf(fid,'\n');
% 
% for k=1:length(t)
%     fprintf(fid,'%.7E',t(k));
%     for j=1:length(real_r_norm)
%         fprintf(fid,' %.7E',u(k,j));        
%     end
%     fprintf(fid,'\n');    
% end