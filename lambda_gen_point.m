function lambda_gen = lambda_gen_point(x,j)
%Extracting the data from the input
N = (length(x)-1)/3;
u = reshape(x(1:3*N), 3, N);
omega = x(end);

u_j = u(:,j);

uj_minus_u = u_j - u;
uj_minus_u(:,j) = [];

del_H = -sum(sum(uj_minus_u.^2,1).^(-1).*uj_minus_u,2) -omega*[0;0;1];

lambda_gen = -dot(del_H,u_j);