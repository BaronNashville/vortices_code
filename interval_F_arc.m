function r_min = interval_F_arc(x, m, poles, u_tilde, approx, tangent, r_star)
%Using first order radii polynomial to rigorously verify solutions to the pseudo-arclength function

x = intval(x);

y = midrad(mid(x),r_star);
Id = eye(length(x));
A = DF_arc(x, m, poles, u_tilde, tangent)^-1;

Y = sup(norm(A*F_arc(x, m, poles, u_tilde, approx, tangent),inf));
Z = sup(norm(Id - A*DF_arc(y, m, poles, u_tilde, tangent),inf));

r_min = Y/(1-Z);
end