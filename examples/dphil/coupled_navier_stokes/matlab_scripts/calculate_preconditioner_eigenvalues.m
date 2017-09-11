% calculate the eigenvalues of the different

system_matrix_test_11

%A = Mat_0x1b86fb0_0;
%A = Mat_0x2805fb0_0;
%A = Mat_0x1a73930_0; % l/d ratio 50
%A = Mat_0x29ea930_0;    % l/d ratio 500
%A = Mat_0x2834930_0;    % l/d ratio 5
%A = Mat_0x26fc020_0;    % asym 1
%A = Mat_0x24b22a0_0;    % asym 2
%A = Mat_0xae92a0_0;     % asym 3
%A = Mat_0x21e1a20_0;    % asym 2x2
A = Mat_0x104f9b0_0;
%A = Mat_0xfa6930_0;

%31 - 27, 4
%size_3d = 2481;
size_3d = 2259; % asym
%size_3d = 4941; % 2x2
%size_3d = 27;
%size_0d = 4;
size_0d = 8;
%size_0d = 16;
size = size_3d + size_0d;
A_3d = A(1:size_3d,1:size_3d);
H = A(size_3d+1:end,1:size_3d);
H_T = A(1:size_3d,size_3d+1:end);
A_0d = A(size_3d+1:end,size_3d+1:end);

id_0d = speye(size_0d);

inv_A_3d = inv(A_3d);
inv_A_0d = inv(A_0d);

P_diagonal = sparse(size,size);
P_diagonal(1:size_3d,1:size_3d) = inv_A_3d;
P_diagonal(size_3d+1:end,size_3d+1:end) = inv_A_0d;


AP_diagonal = A*P_diagonal;


P_schur = sparse(size,size);
P_schur(1:size_3d,1:size_3d) = inv_A_3d;
P_schur(size_3d+1:end,size_3d+1:end) = inv_A_0d;
P_schur(1:size_3d,size_3d+1:end) = -inv_A_3d*H_T*inv_A_0d;


AP_schur = A*P_schur;

%e_diagonal = eigs(AP_diagonal,size)
%e_schur = eigs(AP_schur,size)

e_diagonal = eigs(AP_diagonal,10)
e_schur = eigs(AP_schur,10)


%H*inv(A_3d)*H_T

%full(inv(A_0d))