function [values, err, determinant, tr] = new_stability(x, m, poles)
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

if m >= 3
    if (isintval(x) == 1)
        P = intval(zeros(3*N,1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n));
    else
        P = zeros(3*N,1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n);
    end
    
    for L = 0:floor(m/2)
        if L == 0
            for j = 2:n
                P(:,1+2*(j-2)) = B_hat_fun(j,0);
                %P(:,2+2*(j-2)) = norm(b(1+3*(m-1):3*m,1))^2 *C_hat_fun(j,0) - norm(b(1+3*(m-1):3*m,j))^2*C_hat_fun(1,0)
                P(:,2+2*(j-2)) = norm(w(1))^2 *C_hat_fun(j,0) - norm(w(j))^2*C_hat_fun(1,0);
            end
        elseif L == 1
            P(:,1+2*(n-1)) = z_coord(1)*B_hat_fun(1,1)+1i*C_hat_fun(1,1);
            for j = 2:n
                P(:,2+2*(n-1)+2*(j-2)) = w(j)*B_hat_fun(1,1) - w(1)*B_hat_fun(j,1);
                P(:,3+2*(n-1)+2*(j-2)) = z_coord(j)*w(j)*B_hat_fun(1,1) + 1i*w(1)*C_hat_fun(j,1);
            end
            
            if poles == 1 || poles == -1
                P(:,2+4*(n-1)) = 2*B_hat_fun(1,1)+1i*m*w(1)*(delta_x(:,2)+1i*delta_y(:,2));
            elseif poles == 2
                P(:,2+4*(n-1)) = 2*B_hat_fun(1,1)+1i*m*w(1)*(delta_x(:,2)+1i*delta_y(:,2));
                P(:,3+4*(n-1)) = 2*B_hat_fun(1,1)+1i*m*w(1)*(delta_x(:,1)+1i*delta_y(:,1));
            end
        else
            for j = 1:n
                P(:,2+2*(n-1)+abs(poles) + 2*n*(L-2) + 2*(j-1)) = B_hat_fun(j,L);
                P(:,3+2*(n-1)+abs(poles) + 2*n*(L-2) + 2*(j-1)) = C_hat_fun(j,L);
            end
        end  
    
    end
elseif m == 2
    if (isintval(x) == 1)
        P = intval(zeros(3*N,4*(n-1) + 2*abs(poles) + 2*(floor(m/2)-1) * n));
    else
        P = zeros(3*N,4*(n-1) + 2*abs(poles) + 2*(floor(m/2)-1) * n);
    end
    
    for L = 0:floor(m/2)
        if L == 0
            for j = 2:n
                P(:,1+2*(j-2)) = B_hat_fun(j,0);
                %P(:,2+2*(j-2)) = norm(b(1+3*(m-1):3*m,1))^2 *C_hat_fun(j,0) - norm(b(1+3*(m-1):3*m,j))^2*C_hat_fun(1,0)
                P(:,2+2*(j-2)) = norm(w(1))^2 *C_hat_fun(j,0) - norm(w(j))^2*C_hat_fun(1,0);
            end
        elseif L == 1
            for j = 2:n
                P(:,1+2*(n-1)+2*(j-2)) = z_coord(1)*real(w(1)*conj(w(j)))*B_hat_fun(1,1) -z_coord(1)*norm(w(1))^2 *B_hat_fun(j,1) -imag(w(1)*conj(w(j)))*C_hat_fun(1,1);
                P(:,2+2*(n-1)+2*(j-2)) = z_coord(1)*z_coord(j)*imag(w(1)*conj(w(j)))*B_hat_fun(1,1) - z_coord(1)*norm(w(1))^2 *C_hat_fun(j,1) +z_coord(j)*real(w(1)*conj(w(j)))*C_hat_fun(1,1);
            end
            
            if poles == 1 || poles == -1
                P(:,1+4*(n-1)) = z_coord(1)*imag(w(1))*B_hat_fun(1,1) + real(w(1))*C_hat_fun(1,1) - 2*z_coord(1)*norm(w(1))^2*delta_x(:,2);
                P(:,2+4*(n-1)) = z_coord(1)*real(w(1))*B_hat_fun(1,1) - imag(w(1))*C_hat_fun(1,1) - 2*z_coord(1)*norm(w(1))^2*delta_y(:,2);
            elseif poles == 2
                P(:,1+4*(n-1)) = z_coord(1)*imag(w(1))*B_hat_fun(1,1) + real(w(1))*C_hat_fun(1,1) - 2*z_coord(1)*norm(w(1))^2*delta_x(:,2);
                P(:,2+4*(n-1)) = z_coord(1)*real(w(1))*B_hat_fun(1,1) - imag(w(1))*C_hat_fun(1,1) - 2*z_coord(1)*norm(w(1))^2*delta_y(:,2);
                P(:,3+4*(n-1)) = z_coord(1)*imag(w(1))*B_hat_fun(1,1) + real(w(1))*C_hat_fun(1,1) - 2*z_coord(1)*norm(w(1))^2*delta_x(:,1);
                P(:,4+4*(n-1)) = z_coord(1)*real(w(1))*B_hat_fun(1,1) - imag(w(1))*C_hat_fun(1,1) - 2*z_coord(1)*norm(w(1))^2*delta_y(:,1);
            end
        else
            for j = 1:n
                P(:,1+4*(n-1)+2*abs(poles)+ 2*n*(L-2) + 2*(j-1)) = B_hat_fun(j,L);
                P(:,2+4*(n-1)+2*abs(poles)+ 2*n*(L-2) + 2*(j-1)) = C_hat_fun(j,L);
            end
        end  
    
    end
end

%{
P = [
    u(3,1)*B_hat_fun(1,1)+1i*C_hat_fun(1,1), 2*B_hat_fun(1,1)+1i*3*w(1)*(delta_x(:,2)+1i*delta_y(:,2))
    ];
%}
%{
p1 = z_coord(1)*real(w(1)*conj(w(2)))*B_hat_fun(1,1) -z_coord(1)*norm(w(1))^2 *B_hat_fun(2,1) -imag(w(1)*conj(w(2)))*C_hat_fun(1,1);
p2 = z_coord(1)*z_coord(2)*imag(w(1)*conj(w(2)))*B_hat_fun(1,1) - z_coord(1)*norm(w(1))^2*C_hat_fun(2,1)+z_coord(2)*real(w(1)*conj(w(2)))*C_hat_fun(1,1);


p3 = B_hat_fun(2,0);
p4 = norm(b(4:6,1))^2 *C_hat_fun(2,0) - norm(b(4:6,2))^2*C_hat_fun(1,0);

P = [p1,p2,p3,p4];
%}
    
f_no_sym(new_input);
Q_omega = P' * Df_no_sym(new_input) * P;
values = real(eig(mid(Q_omega)));
determinant = det(mid(Q_omega));
tr = trace(mid(Q_omega));

neg = 0;
zero = 0;
pos = 0;

for k = 1:length(values)
    value = values(k);
    if norm(value) < 1e-10
        zero = zero + 1;
    elseif value > 0
        pos = pos + 1;
    elseif value < 0
        neg = neg + 1;
    end
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


