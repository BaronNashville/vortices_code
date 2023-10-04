function DF = DF_arc(x, m, poles, u_tilde, tangent)
%Computes the Jacobian of the pseudo-arclength map

DF = [tangent';Df(x, m, poles, u_tilde)];
end