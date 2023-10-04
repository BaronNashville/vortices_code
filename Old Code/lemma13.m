%Initialization
clear, clc

%Tetrahedron with vertex as pole
x_rot = 0;
y_rot = 0.955316;
poles = 1;
n = 1;
m = 3;    
coords_1_1_3 = [
    0.9428   -0.4714   -0.0000   -0.4714
         0    0.8165         0   -0.8165
   -0.3333   -0.3333    1.0000   -0.3333
   ];
u = [coords_1_1_3(:,1)];


alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, m, poles);

J_1 = [
    0    0      0
    0    0     -1
    0    1      0
   ];

J_2 = [
    0    0     -1
    0    0      0
    1    0      0
   ];
J_3 = [
    0   -1      0
    1    0      0
    0    0      0
   ];

N = n*m + abs(poles)

X = [u;lambda;alpha;omega]
fixed_X = newton_f(X, m, poles, u_tilde)

X_no_sym = input_no_sym(fixed_X, m, poles)

counter = 3;

X_prime = [-Df_U(fixed_X, m, poles, u_tilde)^-1 * delf_omega(fixed_X, m); 1]

X_prime_no_sym = input_prime_no_sym(X_no_sym, X_prime, m, poles)

while X_no_sym(1)*X_no_sym(5) - X_no_sym(4)*X_no_sym(2)< 1e-5 && counter <= N
    
    tmp = X_no_sym(3:6);
    X_no_sym(3:6) = X_no_sym(1+3*(counter-1):3*counter);
    X_no_sym(1+3*(counter-1):3*counter) = tmp;
    
    tmp = X_prime_no_sym(3:6);
    X_prime_no_sym(3:6) = X_prime_no_sym(1+3*(counter-1):3*counter);
    X_prime_no_sym(1+3*(counter-1):3*counter) = tmp;  
    
    
    
    tmp = X_no_sym(3*N + 2);
    X_no_sym(3*N + 2) = X_no_sym(3*N + counter);
    X_no_sym(3*N + counter) = tmp;
    
    tmp = X_prime_no_sym(3*N + 2);
    X_prime_no_sym(3*N + 2) = X_prime_no_sym(3*N + counter);
    X_prime_no_sym(3*N + counter) = tmp;
    
    counter = counter + 1;
end

lambda_prime_no_sym = X_prime_no_sym(3*N+1:4*N);

a = reshape(X_no_sym(1:3*N), 3, N);
a_prime = reshape(X_prime_no_sym(1:3*N), 3, N);
b = zeros(3,N);
c = zeros(3,N);

B = @(w, w_prime) -2*eye(3)*dot(w, w_prime)/(norm(w)^4) + (8*dot(w,w_prime)/(norm(w))^6)*(w*w') -(2/(norm(w)^4))*(w_prime*w' +w*w_prime');

K = zeros(3*N,3*N);
for row = 1:N
    for col = 1:N
        
        if row == col
            sum = 0;
            for j = 1:N
                if j ~= row
                    sum = sum - B(a(1+3*(row-1):3*row)-a(1+3*(j-1):3*j), a_prime(1+3*(row-1):3*row)-a_prime(1+3*(j-1):3*j));
                end
            end
            K(1+3*(row-1):3*row,1+3*(col-1):3*col) = sum;
        else
            K(1+3*(row-1):3*row,1+3*(col-1):3*col) = B(a(1+3*(row-1):3*row)-a(1+3*(col-1):3*col), a_prime(1+3*(row-1):3*row)-a_prime(1+3*(col-1):3*col));
        end
        
    end
end

zeta_prime = zeros(3*N,3*N);

for i = 1:N
    zeta_prime(1+3*(i-1):3*i,1+3*(i-1):3*i) = lambda_prime_no_sym(i)*eye(3);
end

R_0 = zeros(3*N, 2);

R_0(:,1) = reshape(J_1 * a, 1, 3*N);
R_0(:,1) = reshape(J_2 * a, 1, 3*N);


C = R_0' * (K + zeta_prime) * R_0




