function f = f_eig_pair(vector_value_alpha, k, vk, v_prime, Q)
x = (length(vector_value_alpha)-2);

v = vector_value_alpha(1:x);
value = vector_value_alpha(end-1);
alpha = vector_value_alpha(end);

Id = eye(x);

f = [
    (Q-value*Id)*v+alpha*v
    v(k)-vk
    dot(v, v_prime)+alpha*value
    ];
    