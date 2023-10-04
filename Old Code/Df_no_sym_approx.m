function Df_no_sym_approx = Df_no_sym_approx(x)
%Approximation of the jacobian using finite differencing
k = length(x);
N = (length(x)-1)/4;
Id = eye(k);
h = 1e-8;
Df_no_sym_approx = zeros(3*N, 3*N);
for i = 1:3*N
    Df_no_sym_approx(:, i) = (1/h) * (f_no_sym(x + h*Id(:, i)) - f_no_sym(x));
end