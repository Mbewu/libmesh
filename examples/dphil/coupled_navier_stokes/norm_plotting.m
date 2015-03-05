% do the noem dance

%p2p1
p2p1_data = [0	0.204709	0.526912	0.15936
0	0.0448014	0.211006	0.060766
0	0.00665055	0.0737036	0.0202213
0	0.00112827	0.0294223	0.00972404
0	0.00025998	0.014437	0.00450239];

p1p1_data = [0	0.558678	0.917725	1.79164
0	0.257823	0.621549	1.09716
0	0.074205	0.343281	0.505324
0	0.017617	0.172961	0.237226
0	0.00365308	0.081232	0.0916245];

h = [0.2 0.1 0.05 0.025 0.0125];
h_inv = 1./h;


% figure
% subplot(2,3,1);
% loglog(h_inv,p2p1_data(:,2),'r*');
% hold on
% loglog(h_inv,h.^3,'b');
% loglog(h_inv,h.^2,'g');
% title('p2p1 velocity l2 norm');
% hold off
% 
% subplot(2,3,2);
% loglog(h_inv,p2p1_data(:,3),'r*');
% hold on
% loglog(h_inv,h.^2,'b');
% loglog(h_inv,h.^1,'g');
% title('p2p1 velocity h1 seminorm');
% hold off
% 
% 
% subplot(2,3,3);
% loglog(h_inv,p2p1_data(:,4),'r*');
% hold on
% loglog(h_inv,h.^2,'b');
% loglog(h_inv,h.^1,'g');
% title('p2p1 pressure l2 norm');
% hold off
set(0,'DefaultAxesFontSize',16)

subplot(1,3,1);
loglog(h_inv,p1p1_data(:,2),'b-*');
hold on
loglog(h_inv,h.^2,'r-*');
title('p1p1 velocity l2 norm');
xlabel('1 / h');
ylabel('error');
legend('error','theoretical error');
hold off

subplot(1,3,2);
loglog(h_inv,p1p1_data(:,3),'b-*');
hold on
loglog(h_inv,h,'r-*');
title('p1p1 velocity h1 seminorm');
xlabel('1 / h');
ylabel('error');
legend('error','theoretical error');
hold off


subplot(1,3,3);
loglog(h_inv,p1p1_data(:,4),'b-*');
hold on
loglog(h_inv,h,'r-*');
title('p1p1 pressure l2 norm');
xlabel('1 / h');
ylabel('error');
legend('error','theoretical error');
hold off