function Df_approx = Df_approx(x, m, poles, u_tilde)
%Approximation of the Jacobian of F defined in eq 3.15 using finite differencing
%Used in testing to make sure our analytical Jacobia is correct

k = length(x);
Id = eye(k);
h = 1e-8;
Df_approx = zeros(k-1, k);
for i = 1:k
    Df_approx(:, i) = (1/h) * (f(x + h*Id(:, i), m, poles, u_tilde) - f(x, m, poles, u_tilde));
end