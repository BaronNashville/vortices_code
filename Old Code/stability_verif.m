function [values, bounds, err] = stability_verif(x, m, poles)
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

new_input = input_no_sym(x, m, poles);
a = reshape(new_input(1:3*N), 3, N);

if isintval(a) == 1
    b = intval(zeros(3, N));
elseif isaffari(a) == 1
    b = affari(zeros(3, N));
else
    b = zeros(3, N);
end
if isintval(a) == 1
    c = intval(zeros(3, N));
elseif isaffari(a) == 1
    c = affari(zeros(3, N));
else
    c = zeros(3, N);
end

for j = 1:N
    if norm(a(:,j) - [0;0;1]) < 1e-5 || norm(a(:,j) - [0;0;-1]) < 1e-5
        b(:,j) = [1;0;0];
    else
        b(:,j) = J_3*a(:,j);
    end
    c(:,j) = cross(a(:,j), b(:,j));
end

if isintval(a) == 1
    P = intval(zeros(3*N, 2*N-4));
elseif isaffari(a) == 1
    P = affari(zeros(3*N, 2*N-4));
else
    P = zeros(3*N, 2*N-4);
end

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

Q_omega = intval(P' * Df_no_sym(new_input) * P);

[vectors,values] = eig(mid(Q_omega));
values = real(diag(values));
bounds = zeros(2, length(values));
%{
neg = 0;
ambiguous = 0;
pos = 0;

for i = 1:length(values)
    vector = vectors(:,i);
    value = values(i);
    
    k = 1;
    vk = vector(1);
    
    for j = 1:length(values)
        if vector(j) > vk
            k = j;
            vk = vector(j);
        end
    end
    bounds(:,i) = interval_eig([vector;value], k, vk, Q_omega);
    
    if isnan(bounds(1,i)) || bounds(1,i) < 0
        ambiguous = ambiguous + 1;
    elseif value + bounds(1,i) > 0 && value - bounds(1,i) > 0
        pos = pos + 1;
    elseif value + bounds(1,i) < 0 && value - bounds(1,i) < 0
        neg = neg + 1;
    else
        ambiguous = ambiguous + 1;
    end
end

if ambiguous ~= 0
    err = 3;    %numerical uncertainty - magenta
elseif neg == 0
    err = 0;    %Lyapunov stable - green
elseif mod(neg,2) == 1
    err = 1;    %Unstable solution - red
elseif mod(neg,2) == 0 || mod(pos, 2) == 0
    err = 2;    %inconclusive test - black
end
%}
err = 0;
%{
lambdas1 = zeros(N,1);
lambdas2 = zeros(N,1);
mus1 = zeros(N,1);
mus2 = zeros(N,1);
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
    err = 3;         %Numerical uncertainty - magenta
elseif determinant > 0 && tr > 0
    err = 0;         %Lyapunov stable - green
elseif determinant < 0 && tr < 0
    err = 1;         %unstable - red
else
    err = 2;         %inconclusive - black
end
%}
