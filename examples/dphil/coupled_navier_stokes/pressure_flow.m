%plot some pressure stuff

flow = 60*[0.01 0.05 0.1 0.2 0.3 0.4 0.5 0.75 1.0 1.25 1.5 1.75 2.0];

poiseuille = 0.01019*[0.48716 2.4358 4.8716 9.7432 14.6148 19.4864 24.358 36.537 48.716 60.895 73.074 85.253 97.432];

pedley = 0.01019*[0.809942 5.92805 14.2993 35.1866 60.6368 90.0054 122.845 218.15 329.868 455.989 595.112 746.149 908.215];

ertbruggen = 0.01019*[0.779724 5.72622 13.9229  34.6219 59.9985 89.271 121.982 216.944 328.282 454.009 592.709 743.316 904.935];

reynolds = 0.01019*[0.643898 4.37665 10.582 26.5296 46.351 69.6074 96.0611 175.24 271.884 384.989 513.829 657.835 816.538];

figure
hold on
plot(flow,poiseuille,'r*-');
plot(flow,pedley,'b*-');
plot(flow,ertbruggen,'m*-');
plot(flow,reynolds,'g*-');
legend('poiseuille','pedley','ertbruggen','reynolds');
xlabel('flow rate (L/min)')
ylabel('pressure drop from trachea (cmH2O)')
title('pressure-flow on tree APLE0036266')
