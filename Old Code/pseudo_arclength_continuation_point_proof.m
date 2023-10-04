%Initialization
clear, clc

poles = 1;
n = 1;
m = 3;    
coords_1_1_3 = [
    0.9428   -0.4714   -0.0000   -0.4714
         0    0.8165         0   -0.8165
   -0.3333   -0.3333    1.0000   -0.3333
   ];
u = [coords_1_1_3(:,1)];


alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, m, poles);

steps = 100;
step_size = 1e-4;
max_step_size = 1e-3;
r_star = 1e-14;

solution_proof = 1;
adaptive = 0;
detect_bifurcation = 1;
max_bifurcation_jumps = 100;

initial_solution = [u;lambda;alpha;omega];
%{
fixed_initial_solution = [
    0.1072
    0.9079
    0.4053
   -0.1911
    0.6465
   -0.7386
    0.9852
   -0.0104
    0.1712
    6.3635
    5.3377
    6.1535
   -0.0000
    0.4484
  ];
%}
fixed_initial_solution = newton_f(initial_solution, m, poles, u_tilde);


solutions = zeros(4*n + 2, steps + 1);
solutions(:,1) = fixed_initial_solution;

tangents = zeros(4*n + 2, steps + 1);
nullspace = null(Df(fixed_initial_solution, m, poles, u_tilde));
tangents(:,1) = nullspace(:,1);
if tangents(end,1) < 0
    tangents(:,1) = -tangents(:,1);
end
%{
tangents(:,1) = [
    -0.1123
    0.0407
   -0.0615
    0.2217
   -0.0476
   -0.0991
    0.0343
    0.0235
   -0.1960
   -0.4338
    0.6012
   -0.3357
    0.0000
   -0.4672
   ];
%}
vectorized_solutions = zeros(3, n*m*(steps+1));
vectorized_solutions(:, 1:n*m) = vectorize(fixed_initial_solution, m);

momentum = zeros(1, steps+1);
energy = zeros(1, steps+1);
[momentum(1), energy(1)] = plot_values(solutions(:,1), m, poles);

if poles == 0
    N = n*m;
elseif poles == 1 || poles == -1
    N = n*m+1;
else
    N = n*m+2;
end
values = zeros(2*N-4,steps+1);
[values(:,1), err] = stability(solutions(:,1), m, poles);
if err == 0
    k = [0 1 0];
elseif err == 3
    k = [1 0 1];
elseif err == 2
    k = [0 0 0];
elseif err == 1
    k = [1 0 0];
end

color_delta = 1/(steps);
color1 = [[zeros(n*m,2), ones(n*m,1)];zeros(n*m*steps,3)]; %Coloring continuation plot
color2 = zeros(steps+1,3);                                 %Coloring existence plot
color3 = [k;zeros(steps,3)];                               %Coloring stability plot
r_stars = zeros(1,steps);

if solution_proof == 1
    proof_data = zeros(1,steps);
end

current_bifurcation_jumps = 0;
bifurcation_detected = 0;

normAs = zeros(1, steps);


%Continuation
for i = 1:steps
    %Stepping
    r_stars(i) = r_star;
    disp(['step # = ', num2str(i),', norm of inverse = ' , num2str(norm((DF_arc(solutions(:,i), m, poles,u_tilde, tangents(:,i)))^-1,inf)), ', omega = ', num2str(solutions(end,i))])
    
    approx_solution = solutions(:,i) + step_size * tangents(:,i);    
    [fixed_solution, iter, normAs(i)] = newton_F_arc(approx_solution, m, poles, u_tilde, approx_solution, tangents(:,i));
    
    %Detecting bifurcations
    if detect_bifurcation == 1 && normAs(i) > 2*10^5 %&& i*step_size > 0.1
        [u,s,v] = svd(Df(fixed_solution, m, poles, u_tilde));
        bifurcation_detected = 1;
        current_bifurcation_jumps = current_bifurcation_jumps +1;
    end
    
    vectorized_solutions(:,1 + i*n*m: (i+1)*n*m) = vectorize(fixed_solution, m);
    
    solutions(:,i+1) = fixed_solution;
    nullspace = null(Df(fixed_solution, m, poles, u_tilde));
    tangents(:,i+1) = nullspace(:,1);
    
    %Dealing going along the bifurcation branch
    if bifurcation_detected == 1
        phi_1 = v(:,end);
        phi_2 = v(:,end-1);
        
        if current_bifurcation_jumps < max_bifurcation_jumps
            if dot(phi_2, tangents(:,i+1)) > 1e-2
                new_tangent = phi_1 - dot(phi_1,tangents(:,i+1))/dot(phi_2,tangents(:,i+1)) *phi_2;
            else
                new_tangent = phi_2 - dot(phi_2,tangents(:,i+1))/dot(phi_1,tangents(:,i+1)) *phi_1;
            end
        else
            if norm(dot(phi_1,tangents(:,i+1))) > norm(dot(phi_2, tangents(:,i+1)))
                new_tangent = phi_1;
            else
                new_tangent = phi_2;
            end
        end
        
        tangents(:,i+1) = new_tangent;
        
        %z = fixed_guess; y = new_tangent;
        %save 9_3_bif z y;
        
        bifurcation_detected = 0;
    end
    
    if dot(tangents(:,i+1), tangents(:,i)) < 0
        tangents(:,i+1) = -tangents(:,i+1);
    end
    
    [momentum(i+1), energy(i+1)] = plot_values(solutions(:,i+1), m, poles);
    
    %Verifying solution and stability
    if solution_proof == 1
        proof_data(i) = interval_F_arc(solutions(:,i+1), m, poles, u_tilde, solutions(:,i+1), tangents(:,i+1), r_star);
        if proof_data(i) > 0 && proof_data(i) <= r_star
            solution_with_bounds = infsup(solutions(:,i+1) - proof_data(i),solutions(:,i+1) + proof_data(i));
            [values(:,i+1), tmp, err] = stability_verif(solution_with_bounds, m, poles);
            err = ~stability_prop_point(solutions(:,i+1), m, poles, proof_data(i));
        else
            values(:,i+1) = stability(solutions(:,i+1), m, poles);
            err = 3; %numerical uncertainty
        end
    else
        [values(:,i+1), err] = stability(solutions(:,i+1), m, poles);
    end
    
    %Coloring
    if solution_proof == 0 || (proof_data(i) > 0 && proof_data(i) <= r_star)
        g = [0 0 1] + i*color_delta*[0 1 0];
    else
        g = [0 0 0];
    end
    h = repmat(g, n*m, 1);
    
    if err == 0
        k = [0 1 0];    %green
    elseif err == 1
        k = [1 0 0];    %red
    elseif err == 2
        k = [0 0 0];    %black
    elseif err == 3
        k = [1 0 1];    %magenta
    end
    
    if i == 1
        if solution_proof == 0 || (proof_data(i) > 0 && proof_data(i) <= r_star)
            color1(1:n*m,:) = h - color_delta*[0 1 0];
            color2(1, :) = g - color_delta*[0 1 0];
        else
            color1(1:n*m,:) = h;
            color2(1, :) = g;
        end
    end
    color1(1+n*m*(i):n*m*(i+1),:) = h;
    color2(i+1,:) = g;
    color3(i+1,:) = k;
    
    if solution_proof == 1
        r_star = proof_data(i) + 1e-14;
    end
    
    if adaptive == 1
        disp(step_size)
        step_size = step_size*2^((3-iter)/4);
        step_size = min([step_size max_step_size]);
    end
end


%Graphing

figure(1)
subplot(2,2,1)
[X, Y, Z] = sphere;
surf(X,Y,Z,'FaceColor', [0 0 0], 'EdgeColor', 0.8*[1,1,1], 'FaceAlpha', 0.2);
hold on
scatter3(vectorized_solutions(1,:),vectorized_solutions(2,:),vectorized_solutions(3,:),10, color1,'filled')
if poles == 1
    scatter3(0,0,1, 'k','filled')
elseif poles == -1
    scatter3(0,0,1, 'k', 'filled')
elseif poles == 2
    scatter3(0,0,1, 'k','filled')
    scatter3(0,0,-1, 'k','filled')
end
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

if solution_proof == 1
    subplot(2,2,2)
    scatter(momentum(2:end),proof_data,10,color2(2:end,:),'filled')
    hold on
    scatter(momentum(2:end),r_stars, 2, 'k')    
    set(gca,'FontSize',15)
    axis ([0 inf 0 inf])
    set(gca,'yscale','log')
    xlabel('$$\phi$$', 'Interpreter', 'latex', 'FontSize', 25)
    ylabel('$$r_0$$', 'Interpreter', 'latex', 'FontSize', 25)
    title('Existence Bound Along the Continuation')
end

subplot(2,2,3)
for i = 1:2*N-4
    scatter(momentum, values(i,:), 5, color3, 'filled');
    hold on   
end
plot(momentum, zeros(1,steps+1), 'k')
set(gca,'FontSize',15)
xlabel('$$\phi$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$\lambda$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')

subplot(2,2,4)
scatter(momentum, energy, 5, color3, 'filled');
hold on
set(gca,'FontSize',15)
xlabel('$$\phi$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$H$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')


figure(5)
scatter(solutions(3,:), solutions(6,:), 5, color3, 'filled');
%scatter3(solutions(3,:), solutions(6,:), solutions(9,:), 5, color3, 'filled');
hold on
%plot(-1:0.01:1, -1:0.01:1, 'k');
%plot(-1:0.01:1, 1:-0.01:-1, 'k');
%}

if detect_bifurcation == 1
    figure(6)
    scatter(momentum(2:end), normAs, 5, 'b', 'filled');
    set(gca,'yscale','log')
end
