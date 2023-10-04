function directional = directional(x, x_bar, m, poles, u_tilde)
%Directional derivative in the direction of x_bar of the map F defined in eq 3.15

%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);
lambda = x(3*n + 1: 4*n);
alpha = x(4*n + 1);

u_bar = reshape(x_bar(1:3*n), 3, n);
lambda_bar = x_bar(3*n + 1: 4*n);
alpha_bar = x_bar(4*n + 1);
omega_bar = x_bar(4*n + 2);

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
    directional = intval(zeros(4*n + 1, 1));   
else
    directional = zeros(4*n + 1, 1);
end

%Computing del G_1, ..., G_n
for j = 1:n
    directional(3*n + j, 1) = u(:, j)'*u_bar(:,j);
end

%Computing del sum_{j = 1 to n} {J_3*u_j dot u_tilde_j
directional(end, 1) = sum(diag((J_3 * u_tilde)'*u_bar));

%Computing del f_1, ..., f_n
for j = 1:n
    u_j = u(:, j);
    u_bar_j = u_bar(:,j);

    if isintval(x) == 1
        ring_j = intval(zeros(3,m));    %Contains the position of all vortices in the j'th ring
    else
        ring_j = zeros(3, m);
    end

    for i = 1:m
        ring_j(:,i) = g^i * u_j;
    end

    uj_minus_ring_j = u_j - ring_j;
    uj_minus_ring_j(:, end) = [];

    uj_minus_ring_j_T = uj_minus_ring_j';

    for i = 1:m-1
        uj_minus_ring_j_T(i:end,:) = uj_minus_ring_j_T(i:end,:) * g;
    end
    uj_minus_ring_j_T = uj_minus_ring_j' - uj_minus_ring_j_T;

    sum1 = zeros(3,1);            
    tmp = sum(uj_minus_ring_j.^(2), 1).^(-1);
    for i = 1:m-1
        sum1 = sum1 + tmp(i) * (Id - g^i)*u_bar_j;
    end            
    sum2 = sum(uj_minus_ring_j.^(2), 1).^(-2) .* uj_minus_ring_j*uj_minus_ring_j_T*u_bar_j;


    sum3 = zeros(3,1);
    sum4 = zeros(3,1);


    for j_prime = 1:n
        if j_prime ~= j
            u_j_prime = u(:, j_prime);
            u_bar_j_prime = u_bar(:, j_prime);

            if isintval(x) == 1
                ring_j_prime = intval(zeros(3,m));
                bar_ring_j_prime = intval(zeros(3,m));
            else
                ring_j_prime = zeros(3,m);
                bar_ring_j_prime = zeros(3,m);
            end
            
            for i = 1:m
                ring_j_prime(:,i) = g^i * u_j_prime;
                bar_ring_j_prime(:,i) = g^i * u_bar_j_prime;
            end
            uj_minus_ring_j_prime = u_j - ring_j_prime;
            u_bar_j_minus_bar_ring_j_prime = u_bar_j - bar_ring_j_prime;
            
            sum3 = sum3 + sum(sum(uj_minus_ring_j_prime.^2, 1).^(-1) .* u_bar_j_minus_bar_ring_j_prime, 2);
            sum4 = sum4 + sum(uj_minus_ring_j_prime.^2, 1).^(-2) .* uj_minus_ring_j_prime * sum(uj_minus_ring_j_prime' .* u_bar_j_minus_bar_ring_j_prime', 2);
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
    sum5 = sum5*u_bar_j;
    sum6 = sum6*u_bar_j;

    deljx = -m*(sum1 - 2*sum2 + sum3 - 2*sum4 + sum5 -2*sum6) - m*omega_bar*[0;0;1] + m*lambda_bar(j)*u_j + m*lambda(j)*u_bar_j + alpha_bar*J_3*u_j + alpha*J_3*u_bar_j;

    directional(1 + 3*(j-1):3*j, 1) = deljx;   
end