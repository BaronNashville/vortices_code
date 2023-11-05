function [values, err] = stability_m1(x, m, poles)
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

%Removing all symmetry from the input
new_input = input_no_sym(x, m, poles);

crossed = 0;

%Formatting the data so that the proof works
counter = 3;
while new_input(1)*new_input(5) - new_input(4)*new_input(2)< 1e-2 && counter <= N
    
    tmp = new_input(4:6);
    new_input(4:6) = new_input(1+3*(counter-1):3*counter);
    new_input(1+3*(counter-1):3*counter) = tmp;  
    
    tmp = new_input(3*N + 2);
    new_input(3*N + 2) = new_input(3*N + counter);
    new_input(3*N + counter) = tmp;
    
    counter = counter + 1;
end

a = reshape(new_input(1:3*N), 3, N);

if isintval(x) == 1
    b = intval(zeros(3, N));
else
    b = zeros(3, N);
end
if isintval(x) == 1
    c = intval(zeros(3, N));
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

if isintval(x) == 1
    M = intval(zeros(3*N, 2*N-4));
else
    M = zeros(3*N, 2*N-4);
end

for j = 1:N-2
    %u
    if j == 1
        M(1:6,1) = [cross(a(:,1),cross(b(:,2),c(:,2)));-cross(a(:,1),cross(b(:,2),c(:,2)))];
    else
        M(1:6,j) = [cross(a(:,1),cross(b(:,2),b(:,j+2)));-dot(a(:,1),b(:,j+2))*b(:,2)];
        M(3*(j+2)-2:3*(j+2),j) = dot(a(:,1),b(:,2))*b(:,j+2);
    end
    %v
    M(1:6,N-2+j) = [cross(a(:,1),cross(b(:,2),c(:,j+2)));-dot(a(:,1),c(:,j+2))*b(:,2)];
    M(3*(j+2)-2:3*(j+2),N-2+j) = dot(a(:,1),b(:,2))*c(:,j+2);
end
M_cal = M' * Df_no_sym(new_input) * M;
values = real(eig(mid(M_cal)));

neg = 0;
zero = 0;
pos = 0;

for i = 1:length(values)
    value = values(i);
    if norm(value) < 1e-15
        zero = zero + 1;
    elseif value > 0
        pos = pos + 1;
    elseif value < 0
        neg = neg + 1;
    end
end

for i = 1:length(values)
    for j = 1:length(values)
        vi = values(i);
        vj = values(j);

        if (i ~= j && abs(vi-vj) < 1e-2)
           crossed = 1;
        end

    end
end

if zero ~= 0
    err = 2;    %inconclusive test - black
elseif neg == 0
    err = 0;    %Lyapunov stable - green
elseif mod(neg,2) == 1
    err = 1;    %Unstable solution - red
elseif mod(neg,2) == 0 || mod(pos, 2) == 0
    err = 2;    %inconclusive test - black
end

if crossed
    %err = 4;    %2 eigenvalues have crossed
end

end