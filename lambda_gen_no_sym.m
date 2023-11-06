function lambda_gen = lambda_gen_no_sym(x)
%Computing the values of lambda which will work for our zero finding problem without symmetries
%Extracting the data from the input
N = (length(x)-1)/3;
u = reshape(x(1:3*N), 3, N);
omega = x(end);

e3 = [0;0;1];

if isintval(x) == 1
    lambda_gen = intval(zeros(N,1));
else
    lambda_gen = zeros(N,1);
end

for j = 1:N
    %Extracting u_j
    u_j = u(:,j);

    %uj_minus_u contains the difference between the jth vortex
    %and the position of all vortices
    uj_minus_u = u_j - u;
    
    %we want to sum over all other vortices, so remove the u_j - u_j term
    uj_minus_u(:,j) = [];
    
    %inverse_norm_square_uj_minus_u contains 1 over the 2 norm
    %squared of u_j - u_i in its ith coodrinate    
    inverse_norm_square_uj_minus_u = sum(uj_minus_u.^2,1).^(-1);

    del_H = -sum(inverse_norm_square_uj_minus_u.*uj_minus_u,2);

    lambda_gen(j) = -dot(del_H-omega*e3,u_j);
end