function F = F_arc(x, m, poles, u_tilde, approx, tangent)
%The pseudo-arclength map

F = [dot(x - approx, tangent); f(x, m , poles, u_tilde)];
end