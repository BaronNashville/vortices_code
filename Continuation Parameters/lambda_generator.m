function lambda_gen = lambda_generator(u, omega, m, poles)
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
    
    ring = zeros(3, m);             %Contains the position of all vortices in the jth ring
    for i = 1:m
        ring(1:3,i) = g^i * u_j;
    end
    
    uj_minus_ring = u_j - ring; 
    uj_minus_ring(:,end) = [];
    sum1 = uj_minus_ring*(sum(uj_minus_ring.^(2),1).^(-1))';
    
    sum2 = zeros(3,1);
    
    for j_prime = 1:n
        if j_prime ~= j
            u_j_prime = u(:, j_prime);
            
            ring_prime = zeros(3, m);             %Contains the position of all vortices in the j'th ring
            for i = 1:m
                ring_prime(1:3,i) = g^i * u_j_prime;
            end
    
            u_j_minus_ring_prime = u_j - ring_prime;
            sum2 = sum2 + u_j_minus_ring_prime*(sum(u_j_minus_ring_prime.^(2),1).^(-1))';
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
    
    lambda_gen(j) = dot(sum1 + sum2 + sum3 + omega*u_j, u_j);        
end