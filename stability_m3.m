function [values, err] = stability_m3(x, m, poles)
%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);
lambda = x(3*n + 1: 4*n);
alpha = x(4*n + 1);
omega = x(4*n + 2);
N = n*m + abs(poles);

%Setting the necessary constants
J_3 = [
    0  -1  0
    1   0  0
    0   0  0
    ];

zeta = 2*pi/m;

g = [
    cos(zeta),  sin(zeta), 0;
    -sin(zeta),  cos(zeta),  0;
    0,          0,          1
    ];

a = zeros(3*m, n);
b = zeros(3*m, n);
c = zeros(3*m, n);

new_input = input_no_sym(x, m, poles);

crossed = 0;

for k = 1:m
    a(1+3*(k-1):3*k,:) = expm((k)*J_3*zeta)*u;
    
    b(1+3*(k-1):3*k,:) = J_3 * a(1+3*(k-1):3*k,:);
    
    for j = 1:n
        c(1+3*(k-1):3*k,j) = cross(b(1+3*(k-1):3*k,j),a(1+3*(k-1):3*k,j));
    end
end

B = zeros(3*N*m, n);
C = zeros(3*N*m, n);

for k = 1:m
   for j = 1:n
       B(1+3*N*(k-1)+3*m*(j-1) + 3*(k-1):1+3*N*(k-1)+3*m*(j-1) + 3*(k-1) + 2,j) = b(1+3*(k-1):3*k,j);
       C(1+3*N*(k-1)+3*m*(j-1) + 3*(k-1):1+3*N*(k-1)+3*m*(j-1) + 3*(k-1) + 2,j) = c(1+3*(k-1):3*k,j);
   end
end

B_hat = zeros(3*N*(floor(m/2)+1),n);
C_hat = zeros(3*N*(floor(m/2)+1),n);

for L = 0:floor(m/2)
   for j = 1:n
       B_tmp = reshape(B(:,j), 3*N, m);
       C_tmp = reshape(C(:,j), 3*N, m);
       for k = 1:m
           B_hat(1+3*N*L:3*N*(L+1), j) = B_hat(1+3*N*L:3*N*(L+1), j) + expm(1i*L*k*zeta)*B_tmp(:,k);
           C_hat(1+3*N*L:3*N*(L+1), j) = C_hat(1+3*N*L:3*N*(L+1), j) + expm(1i*L*k*zeta)*C_tmp(:,k);
       end
   end
end

B_hat_fun = @(j,L) B_hat(1+3*N*L:3*N*(L+1), j);
C_hat_fun = @(j,L) C_hat(1+3*N*L:3*N*(L+1), j);

delta_x = zeros(3*N, 2);
delta_x(end-5:end-3,1) = [1;0;0];
delta_x(end-2:end,2) = [1;0;0];

delta_y = zeros(3*N, 2);
delta_y(end-5:end-3,1) = [0;1;0];
delta_y(end-2:end,2) = [0;1;0];

w = @(j) u(1,j) -1i*u(2,j);
x_coord = @(j) u(1,j);
y_coord = @(j) u(2,j);
z_coord = @(j) u(3,j);


if (isintval(x) == 1)
    M = intval(zeros(3*N,1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n));
else
    M = zeros(3*N,1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n);
end

for L = 0:floor(m/2)
    if L == 0
        for j = 2:n
            M(:,1+2*(j-2)) = B_hat_fun(j,0);
            M(:,2+2*(j-2)) = norm(w(1))^2 *C_hat_fun(j,0) - norm(w(j))^2*C_hat_fun(1,0);
        end
    elseif L == 1
        M(:,1+2*(n-1)) = z_coord(1)*B_hat_fun(1,1)+1i*C_hat_fun(1,1);
        for j = 2:n
            M(:,2+2*(n-1)+2*(j-2)) = w(j)*B_hat_fun(1,1) - w(1)*B_hat_fun(j,1);
            M(:,3+2*(n-1)+2*(j-2)) = z_coord(j)*w(j)*B_hat_fun(1,1) + 1i*w(1)*C_hat_fun(j,1);
        end

        if poles == 1 || poles == -1
            M(:,2+4*(n-1)) = 2*B_hat_fun(1,1)+1i*m*w(1)*(delta_x(:,2)+1i*delta_y(:,2));
        elseif poles == 2
            M(:,2+4*(n-1)) = 2*B_hat_fun(1,1)+1i*m*w(1)*(delta_x(:,2)+1i*delta_y(:,2));
            M(:,3+4*(n-1)) = 2*B_hat_fun(1,1)+1i*m*w(1)*(delta_x(:,1)+1i*delta_y(:,1));
        end
    else
        for j = 1:n
            M(:,2+4*(n-1)+abs(poles) + 2*n*(L-2) + 2*(j-1)) = B_hat_fun(j,L);
            M(:,3+4*(n-1)+abs(poles) + 2*n*(L-2) + 2*(j-1)) = C_hat_fun(j,L);
        end
    end  

end

M_cal = M' * Df_no_sym(new_input) * M;

values = zeros(1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n, 1);
weight = zeros(1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n, 1);

neg = 0;
zero = 0;
pos = 0;
for L = 0:floor(m/2)
    if L == 0
        M_cal_L = M_cal(1:2*(n-1),1:2*(n-1));
        current_values = eig(mid(M_cal_L));
        values(1:2*(n-1)) = real(current_values);
        weight(1:2*(n-1)) = ones(2*n-2, 1);
        
        for i = 1:length(current_values)
            for j = 1:length(current_values)
                vi = current_values(i);
                vj = current_values(j);

                if (i ~= j && abs(vi-vj) < 1e-2)
                   crossed = 1;
                end

            end
        end
        
    elseif L == 1
        M_cal_L = M_cal(1+2*(n-1):1 + 4*(n-1) + abs(poles), 1+2*(n-1):1 + 4*(n-1) + abs(poles));
        current_values = eig(mid(M_cal_L));
        values(1+2*(n-1):1 + 4*(n-1) + abs(poles)) = real(current_values);
        weight(1+2*(n-1):1 + 4*(n-1) + abs(poles)) = 2*ones(1 + 2*(n-1) + abs(poles), 1);        
        
        for i = 1:length(current_values)
            for j = 1:length(current_values)
                vi = current_values(i);
                vj = current_values(j);

                if (i ~= j && abs(vi-vj) < 1e-2)
                   crossed = 1;
                end

            end
        end
        
    else
        M_cal_L = M_cal(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1), 2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1));  
        current_values = eig(mid(M_cal_L));
        values(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = real(current_values);
        if L == m/2
            weight(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = ones(2*n, 1);
        else
            weight(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = 2*ones(2*n, 1);            
        end
        
        for i = 1:length(current_values)
            for j = 1:length(current_values)
                vi = current_values(i);
                vj = current_values(j);

                if (i ~= j && abs(vi-vj) < 1e-2)
                   crossed = 1; 
                end

            end
        end
    end
end

for k = 1:length(values)
    value = values(k);
    if norm(value) < 1e-15
        zero = zero + weight(k);
    elseif value > 0
        pos = pos + weight(k);
    elseif value < 0
        neg = neg + weight(k);
    end
end

if abs(x(3) - 0) < 1e-1
    d = 10;
end

if zero ~= 0
    err = 2;    %inconclusive test - magenta
elseif neg == 0
    err = 0;    %Lyapunov stable - green
elseif mod(neg,2) == 1
    err = 1;    %Unstable solution - red
elseif mod(neg,2) == 0 || mod(pos, 2) == 0
    err = 2;    %inconclusive test - magenta
end

if crossed
    %err = 4;    %2 eigenvalues have crossed
end
end