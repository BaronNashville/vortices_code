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
    
    %computing the map by vectorization
    del_H = sum(inverse_norm_square_uj_minus_u.*uj_minus_u,2);
    
    %slotting in our computation in the correct slot
    f(1+(j-1)*3:3*j,1) = -del_H - omega*e3 + lambda(j)*u_j;
end