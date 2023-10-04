n = 5;
m = 7;
N = 3;

x = rand(4*n+2,1);
dir = rand(4*n+2, 1);
u_tilde = rand(3*n,1);
y = rand(4*N+1,1);

poles = 1;

norm(delf_u_readable(x, m, poles, u_tilde) - delf_u(x, m, poles, u_tilde), inf)
%f_readable(x, m, poles, u_tilde)
%Df(x, m, poles, u_tilde)

%norm(f(x, m, poles, u_tilde) - f_readable(x, m, poles, u_tilde), inf)
%norm(Df_approx(x, m, poles, u_tilde) - Df(x, m, poles, u_tilde), inf)
%delf_u(x, m, poles, u_tilde) - delf_u_loops(x, m, poles, u_tilde)

%Df(x, m, poles, u_tilde)*dir - directional(x, dir, m, poles, u_tilde)

%f_no_sym(y)
%Df_no_sym(y) - Df_no_sym_approx(y)


%{
poles = 0;
n = 9;
m = 1;
coords = [
    -0.192975, 0.245317, 0.950042
    -0.098594, 0.0206291, 0.994914
    0.0621589, -0.291795, 0.954459
    -0.175501, -0.243351, 0.953929
    0.108366, -0.0146134, 0.994004
    0.0559616, 0.290386, 0.955272
    -0.32011, -0.0167943, 0.947232
    0.273379, 0.155966, 0.949178
    0.287406, -0.145713, 0.94666
    ]';
%}

%{
poles = 0;
n = 10;
m = 1;
coords = [
    -0.3226, 0.122646, 0.938556
    -0.280504, -0.176247, 0.943533
    -0.112284, 0.301571, 0.946809
    0.330702, 0.046513, 0.942588
    -0.0567831, -0.338791, 0.939146
    0.161413, 0.304671, 0.938681
    0.0531429, -0.150411, 0.987194
    0.0982749, 0.0957603, 0.990541
    -0.131777, 0.0170358, 0.991133
    0.26042, -0.222744, 0.93945
    ]';
%}

%{
poles = 0;
n = 11;
m = 1;
coords = [
    0.176786, 0.0327058, 0.983706
    -0.0308289, -0.455902, 0.889496
    -0.295862, 0.353377, 0.887463
    0.27452, -0.327868, 0.903959
    0.043972, 0.431914, 0.900842
    -0.348784, -0.279335, 0.894607
    -0.432934, 0.0326227, 0.900835
    0.454214, -0.0497866, 0.889501
    -0.114167, 0.136306, 0.984066
    -0.0632152, -0.168303, 0.983706
    0.3363, 0.29427, 0.894599
    ]';
%}

%{
poles = 0;
n = 4;
m = 3;
coords = [
   0.034887632581044   0.136341626351999   0.990047379682701
  -0.249324115175042   0.243756694911171   0.937240715759919
   0.214399606524508   0.302490526779083   0.928726165202127
  -0.042756936922557   0.368292756396256   0.928726165202127
]';
%}

%{
u = reshape(coords, 3*n,1);


alpha = 0;
omega = 50;
u_tilde = u;
lambda = lambda_generator(u, omega, m, poles);
initial_solution = newton_f([u;lambda;alpha;omega], m, poles, u_tilde);

norm(f(initial_solution, m, poles, u_tilde))

sol_bound = interval_f_U(initial_solution, m, poles, u_tilde, 1e-10);

[values, err] = stability_verif(initial_solution, initial_solution, m, poles, sol_bound, 'point');

same_height = 0;

coords = reshape(initial_solution(1:3*n), 3, n);

for i = 1:size(coords,2)
    for j = 1:size(coords,2)
        x = coords(3,i);
        y = coords(3,j);
        
        if i ~= j && abs(x-y) <= 2*sol_bound
            same_height = same_height + 1;
            if (i < j)
                disp(['(', num2str(i), ', ', num2str(j), ')']);
            end
        end
    end
end

vect = vectorize(initial_solution, m);


[X, Y, Z] = sphere;
surf(X,Y,Z,'FaceColor', [0 0 0], 'EdgeColor', 0.8*[1,1,1], 'FaceAlpha', 0.2);
hold on
scatter3(vect(1,:), vect(2,:), vect(3,:), 'filled', 'b')

axis equal
set(gca,'FontSize',15)
xlabel('$$x$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$y$$', 'Interpreter', 'latex', 'FontSize', 25)
zlabel('$$z$$', 'Interpreter', 'latex', 'FontSize', 25)
hXLabel = get(gca,'XLabel');
set(hXLabel,'rotation',0,'VerticalAlignment','middle')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
hZLabel = get(gca,'ZLabel');
set(hZLabel,'rotation',0,'VerticalAlignment','middle')
%}



