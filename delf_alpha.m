function delf_alpha = delf_alpha(x)
%Partial derivative of F defined in eq 3.15 with respect to alpha

%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);

%Setting the necessary constants
J_3 = [
    0,  -1, 0;
    1,  0,  0;
    0,  0,  0
    ];

if isintval(x) == 1
    delf_alpha = intval(zeros(4*n + 1, 1));
else
    delf_alpha = zeros(4*n + 1, 1);
end

%Computing D_alpha G_1, ..., G_n

%Computing D_alpha sum_{j = 1 to n} {J_3*u_j dot u_tilde_j}

%Computing D_alpha f_1, .., f_n
delf_alpha(1:3*n) = reshape(J_3 * u, 3*n, 1);