function lambda_gen = lambda_gen_no_sym(x)
%Extracting the data from the input
N = (length(x)-1)/3;
u = reshape(x(1:3*N), 3, N);
omega = x(end);

if isintval(x) == 1
    lambda_gen = intval(zeros(N,1));
else
    lambda_gen = zeros(N,1);
end
for j = 1:N
    u_j = u(:,j);

    uj_minus_u = u_j - u;
    uj_minus_u(:,j) = [];

    del_H = -sum(sum(uj_minus_u.^2,1).^(-1).*uj_minus_u,2)-omega*[0;0;1];

    lambda_gen(j) = -dot(del_H,u_j);
end