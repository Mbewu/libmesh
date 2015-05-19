clear
% womersley solution
axi = true;
write_file = true
dx = 0.0025;
dt = 0.001;
Re = 1./0.035; %10;
mu = 0.035; %3.5;
P = 1000;
dL = 1;
C = i*P/(dL);       %negative sin solution
C = P/(dL);       %positive cos solution
%C = -P/(dL);       %negative cos solution
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
u0 = C*i/(4*pi) * (1 - cosh(R * sqrt(4*pi/(2*mu)) * (1+i) * r/R) ...
                    /cosh(R * sqrt(4*pi/(2*mu)) * (1+i)));  %n = 1
u0_wom = C*i/(4*pi) * (1 - besselj(0,(R*sqrt(4*pi/mu) * i^(1.5) * r/R)) ...
                    /besselj(0,(R*sqrt(4*pi/mu) * i^(1.5))));       %n = -1
t = 0:dt:1.0;
u = zeros(length(t),length(u0));
u = conj(cos(4*pi*t) + i*sin(4*pi*t))' * u0;   % fucking ' does hermitian transpose
u_wom = conj(cos(4*pi*t) + i*sin(4*pi*t))' * u0_wom;   % fucking ' does hermitian transpose
u = real(u);
u_wom = real(u_wom);

max_mag = max(abs(max(max(real(u)))),abs(min(min(real(u)))));
max_mag_wom = max(abs(max(max(real(u_wom)))),abs(min(min(real(u_wom)))));
max_mag = max(max_mag,max_mag_wom);

%%%%%%%%%%%%%%%%%% plot solution
%u(1,:)
max(abs(u(1,:)))
max(abs(u(2,:)))
max(abs(u(3,:)))

if(~write_file)
    figure
    for j=1:length(t)
       hold off
       plot([-fliplr(real_r) real_r(2:end)],[fliplr(real(u(j,:))) real(u(j,2:end))],'r');
       hold on
       plot([-fliplr(real_r) real_r(2:end)],[fliplr(real(u_wom(j,:))) real(u_wom(j,2:end))],'b');
       legend('parallel plates','womersley')
       axis([-real_r(end) real_r(end) -max_mag max_mag])
       %plot(real(u0(:)));
       title(['t = ' num2str(t(j))])
       pause(0.01)
    end
end

%%%%%%%%%%%%%%%%%% write file
% we want to write it with normalised R
if(write_file)
    filename = 'parallel_plates_solution_axi_cos_negative_times_hundred.dat';
    fid = fopen(filename,'w');
    fprintf(fid,'-1');
    for j=1:length(real_r_norm)
        fprintf(fid,' %.7E',real_r_norm(j));
    end
    fprintf(fid,'\n');

    for k=1:length(t)
        fprintf(fid,'%.7E',t(k));
        for j=1:length(real_r_norm)
            fprintf(fid,' %.7E',u(k,j));        
        end
        fprintf(fid,'\n');    
    end
end