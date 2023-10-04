function Df = Df_eig_pair(vector_value_alpha, k, vk, v_prime, Q)
x = (length(vector_value_alpha)-2);

v = vector_value_alpha(1:x);
value = vector_value_alpha(end-1);
alpha = vector_value_alpha(end);

Id = eye(x);

e_k = zeros(1,x);
e_k(k) = 1;

Df = [
    Q+(alpha-value)*Id                -v                 v
    e_k                       0                  0
    transpose(v_prime)        alpha                  value
    ]; 
    