function Df = Df_no_sym(x)
%Partial derivative with respect to u of our map

%Extracting the data from the input
N = (length(x)-1)/4;
u = reshape(x(1:3*N), 3, N);
lambda = x(3*N + 1: 4*N);
omega = x(end);

Id = eye(3);

if isintval(x) == 1
    Df = intval(zeros(3*N, 3*N));   
else
    Df = zeros(3*N, 3*N);
end

%Computing del f_1, ..., f_n
for j = 1:N
    for y = 1:N
        u_j = u(:, j);
        u_y = u(:, y);
       
        if y == j
            uj_minus_u = u_j - u;
            uj_minus_u(:,j) = [];
            
            sum1 = sum(sum(uj_minus_u.^2).^(-1))*Id;                     
            sum2 = sum(uj_minus_u.^(2), 1).^(-2) .* (uj_minus_u)*uj_minus_u';

            
            deljx = -(sum1 - 2*sum2) + lambda(j)*Id;
        else
            sum1 = -sum((u_j - u_y).^2).^(-1)*Id;
            sum2 = -sum((u_j - u_y).^2).^(-2).*(u_j - u_y)*(u_j - u_y)';
            
            deljx = -(sum1 - 2*sum2);
        end
        
        Df(1 + 3*(j-1):3*j, 1 + 3*(y-1):3*y) = deljx;
    end    
end
    