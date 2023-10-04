function Df = Df(x, m, poles, u_tilde)
%Combines the subfunctions to calculate the Jacobian of F defined in eq 3.15

df_u = delf_u(x, m, poles, u_tilde);
df_lambda = delf_lambda(x, m);
df_alpha = delf_alpha(x);
df_omega = delf_omega(x, m);

Df = [df_u, df_lambda, df_alpha, df_omega];