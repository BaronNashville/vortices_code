function f = f_no_sym(x)
%The unreduced version of eq 3.15

%Extracting the data from the input
N = (length(x)-1)/4;
u = reshape(x(1:3*N), 3, N);
lambda = x(3*N + 1: 4*N);
omega = x(4*N + 1);

e3 = [0;0;1];

if isintval(x) == 1
    f = intval(zeros(3*N, 1));
else
    f = zeros(3*N, 1);
end

%Computing f1, ..., fn
for j = 1:N
    u_j = u(:,j);
    
    uj_minus_u = u_j - u;
    uj_minus_u(:,j) = [];
    
    del_H = sum(sum(uj_minus_u.^2,1).^(-1).*uj_minus_u,2);
    
    f(1+(j-1)*3:3*j,1) = -del_H - omega*e3 + lambda(j)*u_j;
end