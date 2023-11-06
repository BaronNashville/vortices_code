function vectorize = vectorize(x, m)
%Transforms a given input to our main function into something easier to graph

%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);

%Setting the necessary constants
zeta = 2*pi/m;

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

for i = 1:m
    vectorize(1+3*(i-1):3*i,:) = expm(i*J_3*zeta) * u;
end
vectorize = reshape(vectorize,3,n*m);
