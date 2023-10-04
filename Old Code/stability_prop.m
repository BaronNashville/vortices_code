function [prop, determinant, tr] = stability_prop(x_0,x_1, m, poles, bound)
%Extracting the data from the input
s = infsup(0,1);
x = x_0 + s*(x_1 - x_0) + midrad(zeros(length(x_0),1), bound);

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

new_input = input_no_sym(x, m, poles);
a = reshape(new_input(1:3*N), 3, N);

if isintval(x) == 1
    b = intval(zeros(3, N));
    c = intval(zeros(3, N));
    
    P = intval(zeros(3*N, 2*N-4));
else
    b = zeros(3, N);
    c = zeros(3, N);
    
    P = zeros(3*N, 2*N-4);
end
%{
if isintval(x) == 1
    c = intval(zeros(3, N));
else
    c = zeros(3, N);
end
%}

for j = 1:N
    if norm(a(:,j) - [0;0;1]) < 1e-5 || norm(a(:,j) - [0;0;-1]) < 1e-5
        b(:,j) = [1;0;0];
    else
        b(:,j) = J_3*a(:,j);
    end
    c(:,j) = cross(a(:,j), b(:,j));
end

%{
if isintval(x) == 1
    P = intval(zeros(3*N, 2*N-4));
else
    P = zeros(3*N, 2*N-4);
end
%}

for j = 1:N-2
    if j == 1
        P(1:6,1) = [cross(a(:,1),cross(b(:,2),c(:,2)));-cross(a(:,1),cross(b(:,2),c(:,2)))];
    else
        P(1:6,j) = [cross(a(:,1),cross(b(:,2),b(:,j+2)));-dot(a(:,1),b(:,j+2))*b(:,2)];
        P(3*(j+2)-2:3*(j+2),j) = dot(a(:,1),b(:,2))*b(:,j+2);
    end
    
    P(1:6,N-2+j) = [cross(a(:,1),cross(b(:,2),c(:,j+2)));-dot(a(:,1),c(:,j+2))*b(:,2)];
    P(3*(j+2)-2:3*(j+2),N-2+j) = dot(a(:,1),b(:,2))*c(:,j+2);
end

Q_omega = P' * Df_no_sym(new_input) * P;

lambdas1 = intval(zeros(N,1));
lambdas2 = intval(zeros(N,1));
mus1 = intval(zeros(N,1));
mus2 = intval(zeros(N,1));
for j = 2:N
    lambdas1(j) = dot(a(:,3),c(:,j))/dot(b(:,j),b(:,j));
    lambdas2(j) = dot(b(:,3),c(:,j))/dot(b(:,j),b(:,j));
    
    mus1(j) = -dot(a(:,3),b(:,j))/dot(c(:,j),c(:,j));
    mus2(j) = -dot(b(:,3),b(:,j))/dot(c(:,j),c(:,j));
end
v1 = [mus1(2);lambdas1(4:end);mus1(3:end)];
v2 = [mus2(2);lambdas2(4:end);mus2(3:end)];

R = [v1,v2];

S = R' * Q_omega * R;
determinant = S(1,1)*S(2,2) - S(1,2)*S(2,1);
tr = S(1,1) + S(2,2);

if (inf(determinant) < 0 && sup(determinant) > 0) || (inf(tr) < 0 && sup(tr) > 0)
    prop = 3;         %Numerical uncertainty - magenta
elseif determinant > 0 && tr > 0
    prop = 0;         %Lyapunov stable - green
elseif determinant > 0 && tr < 0
    prop = 2;         %inconclusive - black
else
    prop = 1;         %unstable - red
end

determinant = mid(determinant);
tr = mid(tr);
%{
Id = eye(2*N-4); A = mid(Q_omega)^-1;
tmp = (norm(Id - A*Q_omega,inf))
prop = sup(tmp) < 1;
%}