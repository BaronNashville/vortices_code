function input_prime_no_sym = input_prime_no_sym(x_no_sym, x_prime, m, poles)
%Extracting the data from the input
N = (length(x_no_sym)-1)/4;
u_no_sym = reshape(x_no_sym(1:3*N), 3, N);
omega = x_no_sym(4*N + 1);

n = (length(x_prime)-2)/4;
u_prime = reshape(x_prime(1:3*n), 3, n);
lambda_prime = x_prime(3*n + 1: 4*n);
alpha_prime = x_prime(4*n + 1);
omega_prime = x_prime(4*n + 2);

if isintval(x_no_sym) == 1
    lambda_prime_no_sym = intval(zeros(n*m + abs(poles),1)); 
else
    lambda_prime_no_sym = zeros(n*m + abs(poles),1);
end
lambda_prime_no_sym(1:n*m) = repmat(lambda_prime, m, 1);
if poles == 0
    N = n*m;
    u_prime_no_sym = vectorize(x_prime,m);
elseif poles == 1
    N = n*m+1;
    u_prime_no_sym = [vectorize(x_prime,m),[0;0;0]];
    lambda_prime_no_sym(n*m+1) = lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], n*m+1);
elseif poles == -1
    N = n*m;
    u_prime_no_sym = [vectorize(x_prime,m),[0;0;0]];
    lambda_prime_no_sym(n*m+1) = lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], n*m+1);
else
    N = n*m+2;
    u_prime_no_sym = [vectorize(x_prime,m),[0;0;0],[0;0;0]];
    lambda_prime_no_sym(n*m+1) = lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], n*m+1);
    lambda_prime_no_sym(n*m+2) = lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], n*m+2);
end

%{
counter = 3;
while u_no_sym(1,1)*u_no_sym(2,2) - u_no_sym(2,1)*u_no_sym(1,2)< 1e-5
    
    tmp = u_no_sym(:,2);
    u_no_sym(:,2) = u_no_sym(:,counter);
    u_no_sym(:,counter) = tmp;  
    
    tmp = lambda_prime_no_sym(2);
    lambda_prime_no_sym(2) = lambda_prime_no_sym(counter);
    lambda_prime_no_sym(counter) = tmp;
    
    counter = counter + 1;
end
%}

%{
lambda_prime_no_sym
[
lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], 1);
lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], 2);
lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], 3);
lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], 4);
lambda_prime_gen_point([reshape(u_no_sym, 3*N, 1); omega],[reshape(u_prime_no_sym, 3*N, 1); omega_prime], 5)
]
%}
input_prime_no_sym = [reshape(u_prime_no_sym, 3*N,1);lambda_prime_no_sym;omega_prime];