clc
clear all
close all
% Assembling Global Stiffness Matrix from Local Stiffness Matrix
L = input('Enter length of the beam: ');
n_ele = input('Enter the number of elements/members: ');
EI = input('Enter the value of EI (constant for each element): ');
y = input('Enter the numbering of nodes in [1 x n_node] matrix: ');
n_node = n_ele + 1;

% 'x' represents position of each node serially
x = [0 : L/n_ele : L]; % Length of each element is the same
K = zeros((n_node * 2));


for i = 1 : n_ele
    
    p = y(2*i-1:2*i+2);
    Le = x(i + 1) - x(i);
    k_local = EI * [ 12/Le^3   6/Le^2  -12/Le^3   6/Le^2;
                     6/Le^2    4/Le    -6/Le^2    2/Le;
                    -12/Le^3  -6/Le^2   12/Le^3  -6/Le^2;
                     6/Le^2    2/Le    -6/Le^2    4/Le ];
    
    for j = 1 : 4
        for m = 1 : 4
            K(p(j),p(m)) = K(p(j),p(m)) + k_local(j,m);
        end
    end
end
disp(' ');
disp('GLOBAL STIFFNESS MATRIX OF THE BEAM IS: ');
K

n_uk_disp = input('Enter the number of unknown displacements: ');
P_k = input('Enter the known load matrix in the form of [1 x n_uk_disp] matrix: ');
D_k = input('Enter the known displacement matrix in the form of [1 x n_uk_disp] matrix: ');

P_k = P_k';
D_k = D_k';
K_11 = K(1:n_uk_disp,1:n_uk_disp);
K_12 = K(1:n_uk_disp,n_uk_disp+1:n_node*2);
K_21 = K(n_uk_disp+1:n_node*2,1:n_uk_disp);
K_22 = K(n_uk_disp+1:n_node*2,n_uk_disp+1:n_node*2);

D_u = K_11\(P_k - K_12*D_k);
D_u

disp(' ')

P_u = K_21*D_u + K_22*D_k;
P_u

