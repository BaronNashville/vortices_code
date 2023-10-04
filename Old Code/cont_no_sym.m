%Initialization
clear, clc

poles = 1;
n = 5;
m = 2;
coords_1_5_2 = [
    0.4146   -0.4146   -0.9849    0.9849    0.5142   -0.5142    0.4020   -0.4020    0.5188   -0.5188    0
    0.7484   -0.7484   -0.0096    0.0096   -0.8401    0.8401    0.7257   -0.7257   -0.2874    0.2874    0
    0.5176    0.5176    0.1729    0.1729    0.1729    0.1729   -0.5583    0.5583   -0.8051   -0.8051    1
    ];
u = [coords_1_5_2(:,1);coords_1_5_2(:,3);coords_1_5_2(:,7);coords_1_5_2(:,9)];
   

steps = 143;
step_size = 1e-2;
max_step_size = 1e-3;
r_star = 1e-5;

solution_proof = 0;
adaptive = 0;
detect_bifurcation = 1;
max_bifurcation_jumps = 100;

old_input = [
    0.3898
    0.7037
    0.5941
   -0.9560
   -0.0091
    0.2932
    0.4993
   -0.8153
    0.2932
    0.4296
    0.7754
   -0.4628
    0.5834
   -0.3232
   -0.7451
    5.3788
    5.1870
    5.1870
    4.7049
    4.5249
   -0.0000
    0.6376
    ];

[initial_solution, u_tilde] = input_no_sym(old_input, m, poles);
initial_solution = [initial_solution(1:end-1);0;initial_solution(end)];

N = n*m + abs(poles);
m = 1;
poles = 0;

fixed_initial_solution = newton_f(initial_solution, m, poles, u_tilde);

solutions = zeros(4*N + 2, steps + 1);
solutions(:,1) = fixed_initial_solution;

tangents = zeros(4*N + 2, steps + 1);
nullspace = null(Df(fixed_initial_solution, m, poles, u_tilde));
tangents(:,1) = nullspace(:,1);
if tangents(end,1) < 0
    tangents(:,1) = -tangents(:,1);
end

vectorized_solutions = zeros(3, N*m*(steps+1));
vectorized_solutions(:, 1:N*m) = vectorize(fixed_initial_solution, m);

momentum = zeros(1, steps+1);
energy = zeros(1, steps+1);
[momentum(1), energy(1)] = plot_values(solutions(:,1), m, poles);

determinants = zeros(1, steps+1);
traces = zeros(1, steps+1);

values = zeros(2*N-4,steps+1);
[values(:,1), err, determinants(1), traces(1)] = stability(solutions(:,1), m, poles);
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
color1 = [[zeros(N*m,2), ones(N*m,1)];zeros(N*m*steps,3)]; %Coloring continuation plot
color2 = zeros(steps+1,3);                                 %Coloring existence plot
color3 = [k;zeros(steps,3)];                               %Coloring stability plot
r_stars = zeros(1,steps);

if solution_proof == 1
    proof_data = zeros(1,steps);
end

current_bifurcation_jumps = 0;
bifurcation_detected = 0;

normAs = zeros(1, steps);

disp(['step # = ', num2str(0),', norm of inverse = ' , num2str(norm((DF_arc(solutions(:,1), m, poles,u_tilde, tangents(:,1)))^-1,inf)), ', omega = ', num2str(solutions(end,1))])
%Continuation
for i = 1:steps
    %{
    if i == 140
        step_size = 1e-4;
    end
    %}
    %Stepping
    r_stars(i) = r_star;
    
    approx_solution = solutions(:,i) + step_size * tangents(:,i);    
    [fixed_solution, iter, normAs(i)] = newton_F_arc(approx_solution, m, poles, u_tilde, approx_solution, tangents(:,i));
    
    %Detecting bifurcations
    if current_bifurcation_jumps < max_bifurcation_jumps && detect_bifurcation == 1 && normAs(i) > 10^6 %&& i*step_size > 0.1
        [u,s,v] = svd(Df(fixed_solution, m, poles, u_tilde));
        bifurcation_detected = 1;
        current_bifurcation_jumps = current_bifurcation_jumps +1;
    end
    
    vectorized_solutions(:,1 + i*N*m: (i+1)*N*m) = vectorize(fixed_solution, m);
    
    solutions(:,i+1) = fixed_solution;
    nullspace = null(Df(fixed_solution, m, poles, u_tilde));
    tangents(:,i+1) = nullspace(:,1);
    
    %Going along the bifurcation branch
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
        
        load_solution = fixed_solution; load_tangent = new_tangent;
        save p1n5m2_bif_1 load_solution load_tangent;
        
        break
        bifurcation_detected = 0;
        
    end
    
    if dot(tangents(:,i+1), tangents(:,i)) < 0
        tangents(:,i+1) = -tangents(:,i+1);
    end
    
    [momentum(i+1), energy(i+1)] = plot_values(solutions(:,i+1), m, poles);
    
    %Verifying solution and stability
    if solution_proof == 1
        proof_data(i) = interval_pseudo_arclength(solutions(:,i), solutions(:,i+1), m, poles, u_tilde, tangents(:,i), tangents(:,i+1), r_star);
        if proof_data(i) > 0 && proof_data(i) <= r_star
            solution_with_bounds = infsup(solutions(:,i+1) - proof_data(i),solutions(:,i+1) + proof_data(i));
            %solution_with_bounds = solutions(:,i+1);
            [values(:,i+1), tmp, err] = stability_verif(solution_with_bounds, m, poles);
            %{
            if ~stability_prop(solutions(:,i), solutions(:,i+1), m, poles, proof_data(i))
                err = 3;
            end
            %}
            [err, determinants(i+1), traces(i+1)] = stability_prop(solutions(:,i), solutions(:,i+1), m, poles, proof_data(i));
        else
            values(:,i+1) = stability(solutions(:,i+1), m, poles);
            err = 3; %numerical uncertainty
        end
    else
        [values(:,i+1), err, determinants(i+1), traces(i+1)] = stability(solutions(:,i+1), m, poles);
    end
    
    %Coloring
    
    if err == 0
        k = [0 1 0];    %green   - Lyapunov stable
    elseif err == 1
        k = [1 0 0];    %red     - unstable
    elseif err == 2
        k = [0 0 0];    %black   - inconclusive test
    elseif err == 3
        k = [1 0 1];    %magenta - numercial uncertainty
    end
    
    if solution_proof == 0 || (proof_data(i) > 0 && proof_data(i) <= r_star)
        g = [0 0 1] + i*color_delta*[0 1 0];
    else
        g = [0 0 0];
    end
    h = repmat(g, N*m, 1);
    
    if i == 1
        if solution_proof == 0 || (proof_data(i) > 0 && proof_data(i) <= r_star)
            color1(1:N*m,:) = h - color_delta*[0 1 0];
            color2(1, :) = g - color_delta*[0 1 0];
        else
            color1(1:N*m,:) = h;
            color2(1, :) = g;
        end
    end
    color1(1+N*m*(i):N*m*(i+1),:) = h;
    color2(i+1,:) = g;
    color3(i+1,:) = k;
    
    if solution_proof == 1
        r_star = proof_data(i) + step_size/100000;
    end
    
    if adaptive == 1
        disp(step_size)
        step_size = step_size*2^((3-iter)/4);
        step_size = min([step_size max_step_size]);
    end
    
    disp(['step # = ', num2str(i),', norm of inverse = ' , num2str(norm((DF_arc(solutions(:,i+1), m, poles,u_tilde, tangents(:,i+1)))^-1,inf)), ', omega = ', num2str(solutions(end,i+1))])
end


%Graphing

figure(1)
subplot(1,2,1)
[X, Y, Z] = sphere;
surf(X,Y,Z,'FaceColor', [0 0 0], 'EdgeColor', 0.8*[1,1,1], 'FaceAlpha', 0.2);
hold on
scatter3(vectorized_solutions(1,:),vectorized_solutions(2,:),vectorized_solutions(3,:),10, 'b','filled')
%scatter3(vectorized_solutions(1,:),vectorized_solutions(2,:),vectorized_solutions(3,:),10, color1,'filled')
%scatter3(vectorized_solutions(1,1:N),vectorized_solutions(2,1:N),vectorized_solutions(3,1:N),30, color1(1:N,:),'filled')
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
title('Relative Equilibria')

%{
if solution_proof == 1
    subplot(2,2,2)
    scatter(momentum(2:end),proof_data,10,color2(2:end,:),'filled')
    hold on
    scatter(momentum(2:end),r_stars, 2, 'k')    
    set(gca,'FontSize',15)
    axis ([0 inf 0 inf])
    set(gca,'yscale','log')
    xlabel('$$\mu$$', 'Interpreter', 'latex', 'FontSize', 25)
    ylabel('$$r_0$$', 'Interpreter', 'latex', 'FontSize', 25)
    title('Existence Bound Along the Continuation')
end
%}

%{
subplot(2,2,3)
for i = 1:2*N-4
    scatter(momentum, values(i,:), 5, color3, 'filled');
    hold on   
end
plot(momentum, zeros(1,steps+1), 'k')
set(gca,'FontSize',15)
xlabel('$$\mu$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$\lambda$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
%}

subplot(1,2,2)
scatter(momentum, energy, 10, color3, 'filled');
hold on
set(gca,'FontSize',15)
xlabel('$$\mu$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$H$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
title('Energy-Momentum Diagram')


figure(2)
scatter(solutions(3,:), solutions(6,:), 5, color3, 'filled');
%scatter3(solutions(3,:), solutions(6,:), solutions(9,:), 5, color3, 'filled');
hold on
%plot(-1:0.01:1, -1:0.01:1, 'k');
%plot(-1:0.01:1, 1:-0.01:-1, 'k');
%}

if detect_bifurcation == 1
    figure(3)
    scatter(1:steps, normAs, 5, 'b', 'filled');
    set(gca,'yscale','log')
end

%{
figure(4)
subplot(1,2,1)
scatter(momentum, determinants, 5, color3, 'filled');
hold on
set(gca,'FontSize',15)
xlabel('$$\mu$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$det(S)$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')

subplot(1,2,2)
scatter(momentum, traces, 5, color3, 'filled');
hold on
set(gca,'FontSize',15)
xlabel('$$\mu$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$tr(S)$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
%}