function f = f_eig(vector_value, k, vk, M)
%Map G_l defined in eq 5.1

vector = vector_value(1:end-1);
value = vector_value(end);

Id = eye(size(M,1));

f = [(M-value*Id)*vector
    vector(k) - vk];