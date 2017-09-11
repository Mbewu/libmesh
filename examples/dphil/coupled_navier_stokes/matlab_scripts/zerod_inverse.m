% 0D inverse

% equations
% Q_1,in = flux
% Q_1,in - Q_1,out = 0
% P_1,in = P_1,out + R*Q_1,in
% P_1,out = 0

R = 5.;
flux = 1.;
matrix = zeros(4,4);

matrix = [1 0 0 0; 1 -1 0 0; -R 0 1 -1; 0 0 0 1]
           
rhs = [1 0 0 0 ]';

inverse = inv(matrix)

solution = inverse*rhs

%%%
% equations
% Q_1,in - Q_1,out = 0
% P_1,in = P_1,out + R*Q_1,in
% Q_1,in = flux
% P_1,out = 0

matrix = [1 -1 0 0; -R 0 1 -1; 1 0 0 0;  0 0 0 1]
           
rhs = [0 0 1 0 ]';

inverse = inv(matrix)

solution = inverse*rhs
           