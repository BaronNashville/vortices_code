function r_min = interval_pseudo_arclength2(X_0, X_1, m, poles, u_tilde, tangent_0, tangent_1, r_star)
n = (length(X_0)-2)/4;

X_0 = intval(X_0);
X_1 = intval(X_1);
delta_X = X_1 - X_0;

tangent_0 = intval(tangent_0);
tangent_1 = intval(tangent_1);
delta_tangent = tangent_1 - tangent_0;

Id = eye(length(X_0));
s = infsup(0,1);
C = X_0 + s*delta_X + midrad(zeros(length(X_0),1), r_star);
X_s = X_0 + s * delta_X;
tangent_s = tangent_0 + s*delta_tangent;

A = (DF_arc(X_0, m, poles, u_tilde, tangent_0))^-1;
B = A(:, 2:end);

Y1 = sup(norm(B*f(X_0, m, poles, u_tilde), inf));

Y2 = sup(norm(B * Df(X_s, m, poles, u_tilde) * delta_X, inf));

Y = Y1 + Y2;
Z = sup(norm(Id - A*DF_arc(C, m, poles, u_tilde, tangent_s), inf));


r_min = (Y)/(1 - Z);


if r_min <= 0
    r_min = -1;
end
end