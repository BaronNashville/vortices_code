%Initialization

n = length(u)/3;
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, m, poles);

initial_solution = [u;lambda;alpha;omega];
f(initial_solution, m, poles, u_tilde)

steps = 5;
step_size = 0.1;

solutions = zeros(4*n + 2, steps + 1);
solutions(:,1) = initial_solution;

vectorized_solutions = zeros(3, n*m*(steps+1));
vectorized_solutions(:, 1:n*m) = vectorize(initial_solution, m);


%Continuation
for i = 1:steps
    next_omega = solutions(end,i) - step_size;
    approx_solution = solutions(:,i);
    approx_solution(4*n + 2) = next_omega;
    
    fixed_solution = newton_f(approx_solution, m, poles, u_tilde);
    vectorized_solutions(:,1 + i*n*m: (i+1)*n*m) = vectorize(fixed_solution, m);
    
    solutions(:,i+1) = fixed_solution;
end




%Graphing
[X, Y, Z] = sphere;
figure
surf(X,Y,Z,'FaceColor', [0 0 0], 'EdgeColor', 0.8*[1,1,1], 'FaceAlpha', 0.2);
hold on
scatter3(vectorized_solutions(1,:),vectorized_solutions(2,:),vectorized_solutions(3,:), 'b', 'filled')
if poles == 1
    scatter3(0,0,1,'b','filled')
elseif poles == -1
    scatter3(0,0,1, 'b', 'filled')
elseif poles == 2
    scatter3(0,0,1,'b','filled')
    scatter3(0,0,-1,'b','filled')
end
axis equal