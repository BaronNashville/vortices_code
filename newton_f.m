function x = newton_f(x, m, poles, u_tilde)
%Applies the newton algorithm to find a solution of F = 0 defined in 3.15

n = (length(x)-2)/4;
tol = 1e-14;
max_steps = 100;
counter = 0;
while norm(f(x, m, poles, u_tilde)) > tol && counter < max_steps
    x(1:4*n + 1) = x(1:4*n + 1) -(Df_U(x, m, poles, u_tilde)^-1)*f(x, m, poles, u_tilde);
    counter = counter + 1;
end
end