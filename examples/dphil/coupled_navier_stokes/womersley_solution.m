% womersley solution
Re = 10;
P = 800;
R = 0.5;
r = R:-0.001:0;
r = [r fliplr(r(1:end-1))];
u0 = i*P/pi * (1 - besselj(0,(sqrt(pi*Re) * i^(1.5) * r))/besselj(0,(0.5 * sqrt(pi*Re) * i^(1.5))));
t = 0:0.1:10.0;
u = zeros(length(t),length(u0));
u = (cos(pi*t) + i*sin(pi*t))' * u0;
close all
figure
plot(real(u0));

max_mag = max(abs(max(max(real(u)))),abs(min(min(real(u)))));

figure
hold off
for j=1:length(t)
    plot(real(u(j,:)));
    axis([0 length(u0) -max_mag max_mag])
    %plot(real(u0(:)));
    pause(0.5)
end