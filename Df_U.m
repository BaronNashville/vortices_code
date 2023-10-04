function Df_U = Df_U(x, m, poles, u_tilde)
%Combines the subfunctions to calculate the Jacobian of our map
dfu = delf_u(x, m, poles, u_tilde);
dflambda = delf_lambda(x, m);
dfalpha = delf_alpha(x);

Df_U = [dfu, dflambda, dfalpha];