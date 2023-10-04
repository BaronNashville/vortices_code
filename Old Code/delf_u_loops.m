function delf_u = delf_u_loops(x, m, poles, u_tilde)
%Partial derivative with respect to u of our map

%Setting the necessary constants
n = (length(x)-2)/4;
Id = eye(3);
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

%Extracting the data from the input
u = reshape(x(1:3*n), 3, n);
lambda = x(3*n + 1: 4*n);
alpha = x(end-1);

u_tilde = reshape(u_tilde, 3, n); % row j of u_tilde is u_tilde_j

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
       
        if x == j
            
            sum1 = zeros(3,3);            
            for i = 1:m-1
                sum1 = sum1 + (Id - g^i)*norm((Id - g^i)*u_j)^(-2) - 2*(Id - g^i)*u_j*((Id - g^i)*u_j)'*(Id-g^i)*norm((Id - g^i)*u_j)^(-4);
            end
            
            sum2 = zeros(3,3);
            for j_prime = 1:n
                if j_prime ~= j
                    u_j_prime = u(:, j_prime);
                    
                    for i = 1:m
                        sum2 = sum2 + Id*norm(u_j - g^i*u_j_prime)^(-2) - 2*(u_j - g^i*u_j_prime)*(u_j - g^i*u_j_prime)'*norm(u_j - g^i*u_j_prime)^(-4);
                    end
                end
            end
            
            sum5 = zeros(3,3);
            sum6 = zeros(3,3);
            
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
            
            deljx = -m*(sum1 + sum2 + sum5 - 2*sum6) + lambda(j)*Id + alpha*J_3;
        else
            sum1 = zeros(3,3);
            for i = 1:m
                sum1 = sum1 + (-g^i)*norm(u_j - g^i*u_x)^(-2) - 2*(u_j - g^i*u_x)*(u_j - g^i*u_x)'*(-g^i)*norm(u_j - g^i*u_x)^(-4);
            end
            deljx = -m*sum1;
        end
        
        delf_u(1 + 3*(j-1):3*j, 1 + 3*(x-1):3*x) = deljx;
    end    
end
    