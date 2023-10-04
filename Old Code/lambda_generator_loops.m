function lambda_gen = lambda_generator_loops(u, m, poles)
%Generates the starting values of c for the continuation

%Extracting data from the input
n = length(u)/3;
u = reshape(u(1:3*n), 3, n);

%Setting the necessary constants
zeta = 2*pi/m;
g = [
    cos(zeta),  -sin(zeta), 0;
    sin(zeta),  cos(zeta),  0;
    0,          0,          1
    ];

lambda_gen = zeros(n,1);

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
    
    
    lambda_gen(j) = dot(m*(sum1 + sum2 + sum3), u_j);        
end