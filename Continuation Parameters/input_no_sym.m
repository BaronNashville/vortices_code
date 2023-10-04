function [input_no_sym, u_tilde_no_sym] = input_no_sym(x, m, poles)
%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);
lambda = x(3*n + 1: 4*n);
alpha = x(4*n + 1);
omega = x(4*n + 2);

if isintval(x) == 1
    lambda_no_sym = intval(zeros(n*m + abs(poles),1));
 
else
    lambda_no_sym = zeros(n*m + abs(poles),1);
end
lambda_no_sym(1:n*m) = reshape(repmat(lambda', m, 1),1,n*m);
if poles == 0
    N = n*m;
    u_no_sym = vectorize(x,m);
elseif poles == 1
    N = n*m+1;
    u_no_sym = [vectorize(x,m),[0;0;1]];
    lambda_no_sym(n*m+1) = lambda_gen_point([reshape(u_no_sym, 3*N, 1); omega], n*m+1);
elseif poles == -1
    N = n*m;
    u_no_sym = [vectorize(x,m),[0;0;-1]];
    lambda_no_sym(n*m+1) = lambda_gen_point([reshape(u_no_sym, 3*N, 1);omega], n*m+1);
else
    N = n*m+2;
    u_no_sym = [vectorize(x,m),[0;0;1],[0;0;-1]];
    lambda_no_sym(n*m+1) = lambda_gen_point([reshape(u_no_sym, 3*N, 1);omega], n*m+1);
    lambda_no_sym(n*m+2) = lambda_gen_point([reshape(u_no_sym, 3*N, 1);omega], n*m+2);
end

input_no_sym = [reshape(u_no_sym, 3*N,1);lambda_no_sym;omega];
u_tilde_no_sym = u_no_sym;