function r_min = interval_f_U(x, m, poles, u_tilde, r_star)
%Using first order radii polynomial to rigorously verify solutions to the F = 0 defined in 3.15 with omega fixed

y = midrad(x,r_star);
Id = eye(length(x)-1);
A = Df_param(x, m, poles, u_tilde)^-1;

x = intval(x);

Y = sup(norm(A*f(x, m, poles, u_tilde),inf));
Z = sup(norm(Id - A*Df_U(y, m, poles, u_tilde),inf));

r_min = Y/(1-Z);
end