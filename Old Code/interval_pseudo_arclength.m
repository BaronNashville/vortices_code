function r_min = interval_pseudo_arclength(x_0, x_1, m, poles, u_tilde, tangent_0, tangent_1, r_star)
n = (length(x_0)-2)/4;

x_0 = intval(x_0);
x_1 = intval(x_1);
delta_x = x_1 - x_0;

tangent_0 = intval(tangent_0);
tangent_1 = intval(tangent_1);
delta_tangent = tangent_1 - tangent_0;

Id = eye(length(x_0));
s = infsup(0,1);
C = infsup(mid(x_0) - r_star, mid(x_0) + r_star) + s*delta_x;

A = (DF_arc(x_0, m, poles, u_tilde, tangent_0))^-1;
B = A(:, 2:end);

Y1 = sup(norm(B*f(x_0, m, poles, u_tilde), inf));

intFun = @(Tau) directional(x_0 + Tau*s*delta_x, delta_x, m, poles, u_tilde);
steps = 20;
meshValues = intval(zeros(4*n + 1, steps));
for i = 1:steps
    t_int = infsup((i-1)*1/steps, i/steps);
    meshValues(:,i) = 1/steps * intFun(t_int);    
end
Y2 = sup(norm(B*sum(meshValues,2), inf));

Z1 = sup(norm(Id - A*DF_arc(x_0, m, poles, u_tilde, tangent_0), inf));

Z2 = sup(norm(A*[s*delta_tangent';Df(C, m, poles, u_tilde) - Df(x_0, m, poles, u_tilde)], inf));

r_min = (Y1 + Y2)/(1 - Z1 - Z2);


if r_min <= 0
    r_min = 0;
end
end