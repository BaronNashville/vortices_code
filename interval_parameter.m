function r_min = interval_parameter(X_0, X_1, m, poles, u_tilde, r_star)
%Using radii polynomial to rigorously verify a branch of solutions parametrize by omega
%of F = 0 as defined in equation 3.15

X_0 = intval(X_0);
X_1 = intval(X_1);
delta_X = X_1 - X_0;

Id = eye(length(X_0)-1);
s = infsup(0,1);
C = X_0 + s*delta_X + [midrad(zeros(length(X_0)-1,1), r_star); intval(0)];
X_s = X_0 + s * delta_X;

A = Df_param(mid(X_0), m, poles, u_tilde)^-1;

%Computing the bounds Y, Z
Y1 = norm(A*f(X_0, m, poles, u_tilde), inf);
Y2 = norm(A * Df(X_s, m, poles, u_tilde) * delta_X, inf);
Y = Y1 + Y2;

Z = norm(Id - A*Df_param(C, m, poles, u_tilde), inf);

%If Z > 1, return a "bad" interval, else return the minimum radius of existence
if Z > 1
    r_min = -1;
else
    r_min = sup(Y/(1-Z));
end
end