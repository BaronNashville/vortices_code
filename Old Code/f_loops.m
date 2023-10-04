function f = f_loops(x, m, poles, u_tilde)
%Zero finding problem
%Poles = 0 -> no poles
%Poles = 1 -> North pole
%Pole = -1 -> South pole
%Pole = 2 -> Both poles

%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);
lambda = x(3*n + 1: 4*n);
alpha = x(4*n + 1);
omega = x(4*n + 2);
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
e_3 = [0;0;1];

f = zeros(4*n + 1, 1);

%Computing G_1, ..., G_n
G = zeros(n,1);
for j = 1:n
    u_j = u(:, j); 
    G(j) = 1/2 * (norm(u_j)^2 - 1);
end

f(3*n + 1:4*n, 1) = G;

%Computing sum_{j = 1 to n} {J_3*u_j dot u_tilde_j}

transformed_u = J_3 * u;
u_line = reshape(transformed_u, 1, 3*n);
u_tilde_line = reshape(u_tilde, 1, 3*n);

f(end) = sum(u_line .* u_tilde_line);

%Computing f1, ..., fn
for j = 1:n
    u_j = u(:,j);
    
    sum1 = zeros(3,1);
    for i = 1:m-1
        sum1 = sum1 + (u_j - g^i*u_j)/(norm(u_j - g^i*u_j)^2);
    end
    
    sum2 = zeros(3,1);
    
    for j_prime = 1:n
        if j_prime ~= j
            u_j_prime = u(:, j_prime);
            
            for i = 1:m
                sum2 = sum2 + (u_j - g^i*u_j_prime)/(norm(u_j - g^i*u_j_prime)^2);
            end
        end
    end
    
    sum3 = zeros(3,1);
    
    if poles == 1
        sum3 = (u_j - [0;0;1]) * norm(u_j - [0;0;1])^(-2);
    elseif poles == -1
        sum3 = (u_j - [0;0;-1]) * norm(u_j - [0;0;-1])^(-2);
    elseif poles == 2
        north = (u_j - [0;0;1]) * norm(u_j - [0;0;1])^(-2);
        south = (u_j - [0;0;-1]) * norm(u_j - [0;0;-1])^(-2);
        sum3 = north + south;
    end
    
    del_H = -m*(sum1 + sum2 + sum3);        
    f(1 + 3*(j-1): 3*j, 1) = del_H + omega*e_3 + lambda(j)*u_j + alpha*J_3*u_j;
end