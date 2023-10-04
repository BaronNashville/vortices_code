function delf_omega = delf_omega(x, m)
%Partial derivative of F defined in eq 3.15 with respect to omega

%Extracting the data from the input
n = (length(x)-2)/4;

%Setting the necessary constants
e_3 = [0;0;1];

if isintval(x) == 1
    delf_omega = intval(zeros(4*n + 1, 1));
else
    delf_omega = zeros(4*n + 1, 1);
end

%Computing D_omega G_1, ..., G_n

%Computing D_omega sum_{j = 1 to n} {J_3*u_j dot u_tilde_j}

%Computing D_omega f_1, .., f_n
delf_omega(1:3*n) = repmat(-m*e_3,n,1);