function vectorize = vectorize(x, m)
%Transforms a given input to our main function into something easier to graph

%Extracting the data from the input
%n = length(x)/3;
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);

%Setting the necessary constants
zeta = 2*pi/m;
g = [
    cos(zeta),  -sin(zeta), 0;
    sin(zeta),  cos(zeta),  0;
    0,          0,          1
    ];
J_3 = [
    0  -1  0
    1   0  0
    0   0  0
    ];

if isintval(x) == 1
    vectorize = intval(zeros(3*m, n));
else
    vectorize = zeros(3*m, n);
end

%{
if poles == 0
    vectorize = zeros(3, n*m);
elseif poles == 1 || poles == -1
    vectorize = zeros(3, n*m + 1);
else
    vectorize = zeros(3, n*m + 2);
end
%}

for i = 1:m
    vectorize(1+3*(i-1):3*i,:) = expm(i*J_3*zeta) * u;
end
vectorize = reshape(vectorize,3,n*m);

%{
if poles == 1
    vectorize(:, end) = [0;0;1];
elseif poles == -1
    vectorize(:, end) = [0;0;-1];
elseif poles == 2
    vectorize(:, end-1) = [0;0;-1];
    vectorize(:, end) = [0;0;1];
end
%}
