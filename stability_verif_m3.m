function [values, err] = stability_verif_m3(X_0, X_1, m, poles, bound, verif_type)
%Extracting the data from the input

if strcmp(verif_type, 'point')
    X = [midrad(X_0(1:end-1), bound); X_0(end)];
elseif strcmp(verif_type, 'segment') || strcmp(verif_type, 'prop')
    X = hull(X_0,X_1) + [midrad(zeros(length(X_0)-1,1), bound); X_0(end)];
else
    error('Verif type is not valid')
end

n = (length(X)-2)/4;
u = reshape(X(1:3*n), 3, n);
lambda = X(3*n + 1: 4*n);
alpha = X(4*n + 1);
omega = X(4*n + 2);
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

if (isintval(X) == 1)
    a = intval(zeros(3*m, n));
    b = intval(zeros(3*m, n));
    c = intval(zeros(3*m, n));
    
    B = intval(zeros(3*N*m, n));
    C = intval(zeros(3*N*m, n));
    
    B_hat = intval(zeros(3*N*(floor(m/2)+1),n));
    C_hat = intval(zeros(3*N*(floor(m/2)+1),n));
else
    a = zeros(3*m, n);
    b = zeros(3*m, n);
    c = zeros(3*m, n);
    
    B = zeros(3*N*m, n);
    C = zeros(3*N*m, n);
    
    B_hat = zeros(3*N*(floor(m/2)+1),n);
    C_hat = zeros(3*N*(floor(m/2)+1),n);
end

new_input_point = input_no_sym(midrad(X_0, bound), m, poles);
new_input_range = input_no_sym(X, m, poles);

for k = 1:m
    a(1+3*(k-1):3*k,:) = expm((k)*J_3*zeta)*u;
    
    b(1+3*(k-1):3*k,:) = J_3 * a(1+3*(k-1):3*k,:);
    
    for j = 1:n
        c(1+3*(k-1):3*k,j) = cross(b(1+3*(k-1):3*k,j),a(1+3*(k-1):3*k,j));
    end
end

for k = 1:m
   for j = 1:n
       B(1+3*N*(k-1)+3*m*(j-1) + 3*(k-1):1+3*N*(k-1)+3*m*(j-1) + 3*(k-1) + 2,j) = b(1+3*(k-1):3*k,j);
       C(1+3*N*(k-1)+3*m*(j-1) + 3*(k-1):1+3*N*(k-1)+3*m*(j-1) + 3*(k-1) + 2,j) = c(1+3*(k-1):3*k,j);
   end
end

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

if (isintval(X) == 1)
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

M_cal_point = M' * Df_no_sym(new_input_point) * M;
M_cal_range = M' * Df_no_sym(new_input_range) * M;

values = zeros(1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n, 1);
weight = zeros(1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n, 1);
bounds = zeros(2, length(values));
M_cal_Ls = cell(1,1+floor(m/2));

prop = 1;
if strcmp(verif_type, 'point') || strcmp(verif_type, 'segment')
    neg = 0;
    ambiguous = 0;
    pos = 0;

    for L = 0:floor(m/2)
        if L == 0
            M_cal_L = M_cal_point(1:2*(n-1),1:2*(n-1));
            M_cal_Ls{L+1} = M_cal_L;
            [current_vectors, current_values] = eig(mid(M_cal_L));
            current_values = real(diag(current_values));
            values(1:2*n-2) = current_values;
            weight(1:2*n-2) = ones(2*n-2, 1);

            for i = 1:length(current_values)
                vector = current_vectors(:,i);
                value = current_values(i);

                k = 1;
                vk = vector(1);

                for j = 1:length(vector)
                    if abs(vector(j)) > abs(vk)
                        k = j;
                        vk = vector(j);
                    end
                end
                bounds(:,i) = interval_eig([vector;value], k, vk, M_cal_L);            
            end

            M_cal_L_range = M_cal_range(1:2*n-2,1:2*n-2);

            if n > 1 && sum(sum(isnan(M_cal_L_range^-1))) >= 1
                prop = 0;
            end


        elseif L == 1
            M_cal_L = M_cal_point(1+2*(n-1):1 + 4*(n-1) + abs(poles), 1+2*(n-1):1 + 4*(n-1) + abs(poles));
            M_cal_Ls{L+1} = M_cal_L;
            [current_vectors, current_values] = eig(mid(M_cal_L));
            current_values = real(diag(current_values));
            values(1+2*(n-1):1 + 4*(n-1) + abs(poles)) = current_values;
            weight(1+2*(n-1):1 + 4*(n-1) + abs(poles)) = 2*ones(1 + 2*(n-1) + abs(poles), 1);

            for i = 1:length(current_values)
                vector = current_vectors(:,i);
                value = current_values(i);

                k = 1;
                vk = vector(1);

                for j = 1:length(vector)
                    if abs(vector(j)) > abs(vk)
                        k = j;
                        vk = vector(j);
                    end
                end
                bounds(:,2*n-2 + i) = interval_eig([vector;value], k, vk, M_cal_L);            
            end

            M_cal_L_range = M_cal_range(1+2*(n-1):1 + 4*(n-1) + abs(poles), 1+2*(n-1):1 + 4*(n-1) + abs(poles));

            if sum(sum(isnan(M_cal_L_range^-1))) >= 1
                prop = 0;
            end

        else
            M_cal_L = M_cal_point(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1), 2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1));
            M_cal_Ls{L+1} = M_cal_L;
            [current_vectors, current_values] = eig(mid(M_cal_L));
            current_values = real(diag(current_values));
            values(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = current_values;
            if L == floor(m/2)
                weight(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = ones(2*n, 1);
            else
                weight(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = 2*ones(2*n, 1);            
            end

            for i = 1:length(current_values)
                vector = current_vectors(:,i);
                value = current_values(i);

                k = 1;
                vk = vector(1);

                for j = 1:length(vector)
                    if abs(vector(j)) > abs(vk)
                        k = j;
                        vk = vector(j);
                    end
                end
                bounds(:,1+4*(n-1)+abs(poles)+2*n*(L-2)+i) = interval_eig([vector;value], k, vk, M_cal_L);            
            end

            M_cal_L_range = M_cal_range(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1), 2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)); 

            if sum(sum(isnan(M_cal_L_range^-1))) >= 1
                prop = 0;
            end

        end
    end

    for i = 1:length(values)
        if isnan(bounds(1,i)) || bounds(1,i) < 0
            ambiguous = ambiguous + weight(i);
        elseif values(i) + bounds(1,i) > 0 && values(i) - bounds(1,i) > 0
            pos = pos + weight(i);
        elseif values(i) + bounds(1,i) < 0 && values(i) - bounds(1,i) < 0
            neg = neg + weight(i);
        else
            ambiguous = ambiguous + weight(i);
        end
    end

    if prop == 0 || ambiguous ~= 0
        err = 3;    %numerical uncertainty - yellow
    elseif neg == 0
        err = 0;    %Lyapunov stable - green
    elseif mod(neg,2) == 1
        err = 1;    %Unstable solution - red
    elseif mod(neg,2) == 0 || mod(pos, 2) == 0
        err = 2;    %inconclusive test - magenta
    end
    
elseif strcmp(verif_type, 'prop')
    for L = 0:floor(m/2)
        if L == 0
            M_cal_L = M_cal_point(1:2*(n-1),1:2*(n-1));
            M_cal_Ls{L+1} = M_cal_L;
            [current_vectors, current_values] = eig(mid(M_cal_L));
            current_values = real(diag(current_values));
            values(1:2*n-2) = current_values;
            weight(1:2*n-2) = ones(2*n-2, 1);
            
            M_cal_L_range = M_cal_range(1:2*n-2,1:2*n-2);
            if n > 1 && sum(sum(isnan(M_cal_L_range^-1))) >= 1
                prop = 0;
            end
        elseif L == 1
            M_cal_L = M_cal_point(1+2*(n-1):1 + 4*(n-1) + abs(poles), 1+2*(n-1):1 + 4*(n-1) + abs(poles));
            M_cal_Ls{L+1} = M_cal_L;
            [current_vectors, current_values] = eig(mid(M_cal_L));
            current_values = real(diag(current_values));
            values(1+2*(n-1):1 + 4*(n-1) + abs(poles)) = current_values;
            weight(1+2*(n-1):1 + 4*(n-1) + abs(poles)) = 2*ones(1 + 2*(n-1) + abs(poles), 1);

            
            M_cal_L_range = M_cal_range(1+2*(n-1):1 + 4*(n-1) + abs(poles), 1+2*(n-1):1 + 4*(n-1) + abs(poles));
            if sum(sum(isnan(M_cal_L_range^-1))) >= 1
                prop = 0;
            end
        else
            M_cal_L = M_cal_point(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1), 2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1));
            M_cal_Ls{L+1} = M_cal_L;
            [current_vectors, current_values] = eig(mid(M_cal_L));
            current_values = real(diag(current_values));
            values(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = current_values;
            if L == floor(m/2)
                weight(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = ones(2*n, 1);
            else
                weight(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)) = 2*ones(2*n, 1);            
            end

            
            M_cal_L_range = M_cal_range(2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1), 2+4*(n-1)+abs(poles)+2*n*(L-2):1+4*(n-1)+abs(poles)+2*n*(L-1)); 
            if sum(sum(isnan(M_cal_L_range^-1))) >= 1
                prop = 0;
            end
        end
    end
    
    if prop
        err = 0;
    else
        err = 3;
    end
end
end



















