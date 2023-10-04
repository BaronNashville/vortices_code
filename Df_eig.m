function Df = Df_eig(vector_value, k, vk, M)
%Computes the Jacobian of the map G_l defined in eq 5.1

%Extracting the eigenvector and eigenvalue from the input vector_value
vector = vector_value(1:end-1);
value = vector_value(end);

Id = eye(size(M,1));
e_k = zeros(1,size(M,1));
e_k(k) = 1;

Df = [
    M-value*Id   -vector
    e_k            0
    ];