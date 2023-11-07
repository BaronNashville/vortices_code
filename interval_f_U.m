function r_min = interval_f_U(x, m, poles, u_tilde, r_star)
%Using first order radii polynomial to rigorously verify solutions to the F = 0 defined in 3.15 with omega fixed

x = intval(x);

C = midrad(mid(x),r_star);
Id = eye(length(x)-1);
A = Df_param(mid(x), m, poles, u_tilde)^-1;

%Computing the bounds Y, Z
Y = norm(A*f(x, m, poles, u_tilde),inf);
Z = norm(Id - A*Df_U(C, m, poles, u_tilde),inf);

%If Z > 1, return a "bad" interval, else return the minimum radius of existence
if Z > 1
    r_min = -1;
else
    r_min = sup(Y/(1-Z));
end
end