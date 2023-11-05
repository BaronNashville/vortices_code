function r_min = interval_F_arc(x, m, poles, u_tilde, approx, tangent, r_star)
%Using first order radii polynomial to rigorously verify solutions to the pseudo-arclength function

x = intval(x);

C = midrad(mid(x),r_star);
Id = eye(length(x));
A = DF_arc(mid(x), m, poles, u_tilde, tangent)^-1;

%Computing the bounds Y, Z
Y = sup(norm(A*F_arc(x, m, poles, u_tilde, approx, tangent),inf));
Z = sup(norm(Id - A*DF_arc(C, m, poles, u_tilde, tangent),inf));

%If Z > 1, return a "bad" interval, else return the minimum radius of existence
if Z > 1
    r_min = -1;
else
    r_min = Y/(1-Z);
end
end