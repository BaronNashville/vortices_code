function [values, err] = stability_verif_m1(X_0, X_1, m, poles, bound, verif_type)
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

new_input_point = input_no_sym(midrad(X_0, bound), m, poles);
new_input_range = input_no_sym(X, m, poles);

counter = 3;
while new_input_point(1)*new_input_point(5) - new_input_point(4)*new_input_point(2)< 1e-2 && counter <= N
    
    tmp = new_input_point(4:6);
    new_input_point(4:6) = new_input_point(1+3*(counter-1):3*counter);
    new_input_point(1+3*(counter-1):3*counter) = tmp;  
    
    tmp = new_input_point(3*N + 2);
    new_input_point(3*N + 2) = new_input_point(3*N + counter);
    new_input_point(3*N + counter) = tmp;
    
    tmp = new_input_range(4:6);
    new_input_range(4:6) = new_input_range(1+3*(counter-1):3*counter);
    new_input_range(1+3*(counter-1):3*counter) = tmp;  
    
    tmp = new_input_range(3*N + 2);
    new_input_range(3*N + 2) = new_input_range(3*N + counter);
    new_input_range(3*N + counter) = tmp;
    
    counter = counter + 1;
end

a = reshape(new_input_range(1:3*N), 3, N);

if isintval(X) == 1
    b = intval(zeros(3, N));
else
    b = zeros(3, N);
end
if isintval(X) == 1
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

if isintval(X) == 1
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

M_cal_point = M' * Df_no_sym(new_input_point) * M;
M_cal_range = M' * Df_no_sym(new_input_range) * M;

[vectors, values] = eig(mid(M_cal_point));
vectors = real(vectors);
values = real(diag(values));
bounds = zeros(2, length(values));

if strcmp(verif_type, 'point') || strcmp(verif_type, 'segment')
    ambiguous = 0;
    neg = 0;
    pos = 0;
    
    if strcmp(verif_type, 'point')
        d = 10;
    end

    for i = 1:length(values)
        vector = vectors(:,i);
        value = values(i);

        k = 1;
        vk = vector(1);

        for j = 1:length(values)
            if abs(vector(j)) > abs(vk)
                k = j;
                vk = vector(j);
            end
        end
        bounds(:,i) = interval_eig([vector;value], k, vk, M_cal_point);

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

    prop = sum(sum(isnan(M_cal_range^-1))) == 0;

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
    prop = sum(sum(isnan(M_cal_range^-1))) == 0;
    if prop
        err = 0;
    else
        err = 3;
    end
end
end


