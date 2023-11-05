n = 5;
m = 7;
N = 3;

x = rand(4*n+2,1);
dir = rand(4*n+2, 1);
u_tilde = rand(3*n,1);
y = rand(4*N+1,1);

poles = 1;

%norm(delf_u_readable(x, m, poles, u_tilde) - delf_u(x, m, poles, u_tilde), inf)
norm(f_no_sym_readable(x(1:end-1))-f_no_sym(x(1:end-1)))
%Df(x, m, poles, u_tilde)

%norm(f(x, m, poles, u_tilde) - f_readable(x, m, poles, u_tilde), inf)
%norm(Df_approx(x, m, poles, u_tilde) - Df(x, m, poles, u_tilde), inf)
%delf_u(x, m, poles, u_tilde) - delf_u_loops(x, m, poles, u_tilde)

%Df(x, m, poles, u_tilde)*dir - directional(x, dir, m, poles, u_tilde)

%f_no_sym(y)
%Df_no_sym(y) - Df_no_sym_approx(y)




