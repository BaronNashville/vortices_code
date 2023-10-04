function delf_u = delf_u_old(x, m, poles, u_tilde)
%Partial derivative with respect to u of our map

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

delf_u = zeros(4*n + 1, 3*n);

%Computing del G_1, ..., G_n
for j = 1:n
    delf_u(3*n + j, 1 + 3*(j-1): 3*j) = u(:, j)';
end

%Computing del sum_{j = 1 to n} {J_3*u_j dot u_tilde_j
delf_u(end, :) = reshape(J_3' * u_tilde, 1, 3*n);

%Computing del f_1, ..., f_n
for j = 1:n
    for x = 1:n
        u_j = u(:, j);
        u_x = u(:, x);
        ring = zeros(3, m);             %Contains the position of all vortices in the j'th ring
        for i = 1:m
            ring(1:3,i) = g^i * u_j;
        end
       
        if x == j
            ring_minus_uj = ring - u_j;
            ring_minus_uj(:, end) = [];
            
            ring_minus_uj_T = ring_minus_uj';
            
            for i = 1:m-1
                ring_minus_uj_T(i:end,:) = ring_minus_uj_T(i:end,:) * g;
            end
            ring_minus_uj_T = ring_minus_uj_T - ring_minus_uj';
            
            sum1 = zeros(3,3);            
            tmp = sum(ring_minus_uj.^(2), 1).^(-1);
            for i = 1:m-1
                sum1 = sum1 + tmp(i) * (g^i - Id);
            end            
            sum2 = sum(ring_minus_uj.^(2), 1).^(-2) .* (ring_minus_uj)*ring_minus_uj_T;
            
            
            sum3 = zeros(3,3);
            sum4 = zeros(3,3);
            
            
            for j_prime = 1:n
                if j_prime ~=j
                    u_j_prime = u(:, j_prime);
                    ring_minus_uj_prime = ring - u_j_prime;
                    
                    ring_minus_uj_prime_T = ring_minus_uj_prime';
                    for i = 1:m
                        ring_minus_uj_prime_T(i:end,:) = ring_minus_uj_prime_T(i:end,:) * g;
                    end
                    
                    tmp = sum(ring_minus_uj_prime.^2, 1).^(-1);
                    for i = 1:m
                        sum3 = sum3 + tmp(i) * g^i;
                    end
                    sum4 = sum4 + sum(ring_minus_uj_prime.^2, 1).^(-2).*(ring_minus_uj_prime)*(ring_minus_uj_prime_T);
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
            
            deljx = -m*(1/2 * (sum1 -2*sum2) + sum3 - 2*sum4 + sum5 -2*sum6) + lambda(j)*Id + alpha*J_3;
        else
            ring_minus_ux = ring - u_x;
            
            sum1 = sum(sum(ring_minus_ux.^(2), 1).^(-1)) * Id;
            sum2 = sum(ring_minus_ux.^(2), 1).^(-2).*(ring_minus_ux)*(ring_minus_ux)';
            
            deljx = -m*(-sum1 +2*sum2);
        end
        
        delf_u(1 + 3*(j-1):3*j, 1 + 3*(x-1):3*x) = deljx;
    end    
end
    