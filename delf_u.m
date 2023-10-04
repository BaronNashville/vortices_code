function delf_u = delf_u(x, m, poles, u_tilde)
%Partial derivative of F defined in eq 3.15 with respect to u

%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);
lambda = x(3*n + 1: 4*n);
alpha = x(4*n + 1);
u_tilde = reshape(u_tilde, 3, n);


%Setting the necessary constants
zeta = 2*pi/m;
g = [
    cos(zeta),  -sin(zeta), 0;
    sin(zeta),  cos(zeta),  0;
    0,          0,          1
    ];
J_3 = [
    0,  -1, 0;
    1,  0,  0;
    0,  0,  0
    ];
Id = eye(3);

if isintval(x) == 1
    delf_u = intval(zeros(4*n + 1, 3*n));   
else
    delf_u = zeros(4*n + 1, 3*n);
end

%Computing del G_1, ..., G_n
for j = 1:n
    delf_u(3*n + j, 1 + 3*(j-1): 3*j) = u(:, j)';
end

%Computing del sum_{j = 1 to n} {J_3*u_tilde_j dot u_j}
delf_u(end, :) = reshape((J_3 * u_tilde), 1, 3*n);

%Computing del f_1, ..., f_n
for j = 1:n
    for y = 1:n
        u_j = u(:, j);
        u_y = u(:, y);
        
        %Set ring_j to the position of all vortices in the j'th ring
        if isintval(x) == 1
            ring_j = intval(zeros(3,m));
        else
            ring_j = zeros(3, m);
        end
        
        for i = 1:m
            ring_j(:,i) = g^i * u_j;
        end
       
        if y == j
            uj_minus_ring_j = u_j - ring_j;
            uj_minus_ring_j(:, end) = [];
            
            uj_minus_ring_j_T = uj_minus_ring_j';
            
            for i = 1:m-1
                uj_minus_ring_j_T(i:end,:) = uj_minus_ring_j_T(i:end,:) * g;
            end
            uj_minus_ring_j_T = uj_minus_ring_j' - uj_minus_ring_j_T;
            
            sum1 = zeros(3,3);            
            tmp = sum(uj_minus_ring_j.^(2), 1).^(-1);
            for i = 1:m-1
                sum1 = sum1 + tmp(i) * (Id - g^i);
            end            
            sum2 = sum(uj_minus_ring_j.^(2), 1).^(-2) .* (uj_minus_ring_j)*uj_minus_ring_j_T;
            
            
            sum3 = zeros(3,3);
            sum4 = zeros(3,3);
            
            
            for j_prime = 1:n
                if j_prime ~=j
                    u_j_prime = u(:, j_prime);
                    
                    if isintval(x) == 1
                        ring_j_prime = intval(zeros(3,m));   
                    else
                        ring_j_prime = zeros(3,m);
                    end
                    
                    for i = 1:m
                        ring_j_prime(:,i) = g^i * u_j_prime;
                    end
                    uj_minus_ring_j_prime = u_j - ring_j_prime;
                    
                    sum3 = sum3 + sum(sum(uj_minus_ring_j_prime.^2, 1).^(-1)) * Id;
                    sum4 = sum4 + sum(uj_minus_ring_j_prime.^2, 1).^(-2).*(uj_minus_ring_j_prime)*(uj_minus_ring_j_prime)';
                end
            end
            
            
            
            sum5 = 0;
            sum6 = 0;
            
            if poles == 1
                sum5 = norm(u_j - [0;0;1])^(-2) * Id;
                sum6 = (u_j - [0;0;1]) * (u_j - [0;0;1])' * norm(u_j - [0;0;1])^(-4);
            elseif poles == -1
                sum5 = norm(u_j - [0;0;-1])^(-2) * Id;
                sum6 = (u_j - [0;0;-1]) * (u_j - [0;0;-1])' * norm(u_j - [0;0;-1])^(-4);
            elseif poles == 2
                north5 = norm(u_j - [0;0;1])^(-2) * Id;
                north6 = (u_j - [0;0;1]) * (u_j - [0;0;1])' * norm(u_j - [0;0;1])^(-4);
                
                south5 = norm(u_j - [0;0;-1])^(-2) * Id;
                south6 = (u_j - [0;0;-1]) * (u_j - [0;0;-1])' * norm(u_j - [0;0;-1])^(-4);
                
                sum5 = north5 + south5;
                sum6 = north6 + south6;
            end
            
            deljx = -m*(sum1 - 2*sum2 + sum3 - 2*sum4 + sum5 -2*sum6) + m*lambda(j)*Id + alpha*J_3;
        else
            if isintval(x) == 1
                ring_x = intval(zeros(3,m));   
            else
                ring_x = zeros(3,m);
            end
            
            for i = 1:m
                ring_x(:,i) = g^i * u_y;
            end
            uj_minus_ring_x = u_j - ring_x;
            
            uj_minus_ring_x_T = uj_minus_ring_x';            
            for i = 1:m
                uj_minus_ring_x_T(i:end,:) = uj_minus_ring_x_T(i:end,:) * g;
            end
            uj_minus_ring_x_T = - uj_minus_ring_x_T;

            sum1 = zeros(3,3);
            tmp = sum(uj_minus_ring_x.^2, 1).^(-1);
            for i = 1:m
                sum1 = sum1 + tmp(i) * -g^i;
            end
            sum2 = sum(uj_minus_ring_x.^2, 1).^(-2).*(uj_minus_ring_x)*(uj_minus_ring_x_T);
            
            deljx = -m*(sum1 - 2*sum2);
        end
        
        delf_u(1 + 3*(j-1):3*j, 1 + 3*(y-1):3*y) = deljx;
    end    
end
    