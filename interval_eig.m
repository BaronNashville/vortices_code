function int = interval_eig(x, k, vk, M)
%Using radii polynomial to rigorously verify solutions to G_l = 0 defined in eq 5.1

x = intval(x);
M = intval(M);

A = Df_eig(mid(x), k, vk, mid(M))^-1;
Id = eye(length(x));

Y_0 = sup(norm(A*f_eig(x, k, vk, M), inf));
Z_0 = sup(norm(Id - A*Df_eig(x, k, vk, M), inf));
Z_2 = 2*sup(norm(A, inf));

delta = ((1-Z_0)^2 - 4*Z_2*Y_0)^(1/2);

r_min = (1-Z_0 - delta)/(2*Z_2);
r_max = (1-Z_0 + delta)/(2*Z_2);

int = [r_min;r_max];

if imag(delta) ~= 0
    int = [-1; -1];
end
end