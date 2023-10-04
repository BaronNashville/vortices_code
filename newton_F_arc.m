function [x, counter, normA] = newton_F_arc(x, m, poles, u_tilde, approx,tangent)
n = (length(x)-2)/4;
tol = 1e-14;
max_steps = 100;
counter = 0;
while norm(F_arc(x, m, poles, u_tilde, approx,tangent)) > tol && counter < max_steps
    x(1:4*n + 2) = x(1:4*n+2) -(DF_arc(x, m, poles, u_tilde, tangent)^-1)*F_arc(x, m, poles, u_tilde, approx,tangent);
    counter = counter + 1;
    %norm(F_arc(x, m, poles, u_tilde, approx,tangent))
end
normA = norm(DF_arc(x, m, poles, u_tilde, tangent)^-1,inf);
end