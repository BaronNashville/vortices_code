function delf_lambda = delf_lambda(x, m)
%Partial derivative of F defined in eq 3.15 with respect to lambda

%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);

if isintval(x) == 1
    delf_lambda = intval(zeros(4*n + 1, n));
else
    delf_lambda = zeros(4*n + 1, n);
end

%Computing D_lambda G_1, ..., G_n

%Computing D_lambda sum_{j = 1 to n} {J_3*u_j dot u_tilde_j}

%Computing D_lambda f_1, ..., f_n
for j = 1:n
    delf_lambda(1 + 3*(j-1): 3*j, j) = m*u(:, j);
end