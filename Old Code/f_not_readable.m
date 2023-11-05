function f = f(x, m, poles, u_tilde)
%The map F defined in eq 3.15

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

if isintval(x) == 1
    f = intval(zeros(4*n + 1, 1));
else
    f = zeros(4*n + 1, 1);
end

%Computing G_1, ..., G_n
if isintval(x) == 1
    G = intval(zeros(n,1));
else
    G = zeros(n,1);
end
for j = 1:n
    u_j = u(:, j); 
    G(j) = 1/2 * (norm(u_j)^2 - 1);
end

f(3*n + 1:4*n, 1) = G;

%Computing sum_{j = 1 to n} {J_3*u_tilde_j dot u_j}

transformed_u_tilde = J_3 * u_tilde;
u_line = reshape(u, 1, 3*n);
u_tilde_line = reshape(transformed_u_tilde, 1, 3*n);

f(end) = sum(u_line .* u_tilde_line);

%Computing f1, ..., fn
for j = 1:n
    u_j = u(:,j);
    
    if isintval(x) == 1
        ring = intval(zeros(3, m));             %Contains the position of all vortices in the jth ring
    else
        ring = zeros(3, m);
    end
    
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
            
            if isintval(x) == 1
                ring_prime = intval(zeros(3,m));
            else
                ring_prime = zeros(3, m);
            end
            
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
    
    del_H = -m*(sum1 + sum2 + sum3);        
    f(1 + 3*(j-1): 3*j, 1) = del_H - m*omega*e_3 + m*lambda(j)*u_j + alpha*J_3*u_j;
end