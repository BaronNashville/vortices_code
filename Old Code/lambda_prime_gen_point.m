function lambda_gen = lambda_prime_gen_point(x,x_prime,j)
%Extracting the data from the input
N = (length(x)-1)/3;
u = reshape(x(1:3*N), 3, N);
omega = x(end);

u_prime = reshape(x_prime(1:3*N), 3, N);
omega_prime = x_prime(end);

u_j = u(:,j);
u_j_prime = u_prime(:,j);

uj_minus_u = u_j - u;
uj_minus_u(:,j) = [];

uj_prime_minus_u_prime = u_j_prime - u_prime;
uj_prime_minus_u_prime(:,j) = [];

del_H = -sum(sum(uj_minus_u.^2,1).^(-1).*uj_minus_u,2) -omega*[0;0;1];

del_H_prime = -sum(sum(uj_minus_u.^2,1).^(-1).*uj_prime_minus_u_prime,2) + 2*sum(sum(uj_minus_u.^2,1).^(-2).*uj_minus_u*uj_minus_u'*uj_prime_minus_u_prime,2) - [0;0;1];

lambda_gen = (-dot(del_H_prime,u_j) - dot(del_H,u_j_prime));