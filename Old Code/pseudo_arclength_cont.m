%Initialization

starting_step_size = 1e-3;
step_size = starting_step_size;
max_step_size = 1e-3;
r_star = 1e-5;

detect_bifurcation = 1;
max_bifurcation_jumps = 0;

if existence_proof == 1
    proof_data_uniform = zeros(1,steps);
    proof_data_point = zeros(1,steps+1);
end

fixed_initial_solution = newton_f(initial_solution, m, poles, u_tilde);
if existence_proof
    proof_data_point(1) = interval_f_U(fixed_initial_solution, m, poles, u_tilde, 1e-10);
end

solutions = zeros(4*n + 2, steps + 1);
solutions(:,1) = fixed_initial_solution;

tangents = zeros(4*n + 2, steps + 1);

if norm(load_tangent) < 1e-2
    nullspace = null(Df(fixed_initial_solution, m, poles, u_tilde));
    tangents(:,1) = nullspace(:,1);
    if tangents(end,1) < 0
        tangents(:,1) = -tangents(:,1);
    end
else
    tangents(:,1) = load_tangent;
end

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

if m == 1
    values = zeros(2*N - 4, steps+1);
elseif m == 2
    values = zeros(2*(floor(m/2)+1)*(n-1)+2*abs(poles), steps+1);
    
elseif m >=3
    values = zeros(1 + 4*(n-1) + abs(poles) + 2*(floor(m/2)-1) * n, steps+1);
end

if stability_proof
    if (0 <= proof_data_point(1) && proof_data_point(1) <= r_star)
        [values(:,1), err] = stability_verif(solutions(:,1), solutions(:,1), m, poles, proof_data_point(1), 'point');
    else
        values(:,1) = stability(solutions(:,1), m, poles);
        err = 3; %numerical uncertainty
    end
else
    [values(:,1), err] = stability(solutions(:,1), m, poles);
end

k = error_code_to_colour(err);

if segment == 1
    segment_color = [0 0 1];
elseif segment == 2
    segment_color = [1 0 1];
elseif segment == 3
    segment_color = [1, 0.5, 0.5];
elseif segment == 4
    segment_color = [0 1 1];
elseif segment == 5
    segment_color = [0.4940 0.1840 0.5560];
end


color_delta = 1/(steps);
continuation_color = [[repmat(segment_color, n*m, 1)];zeros(n*m*steps,3)]; %Coloring continuation plot
existence_color = zeros(steps+1,3);                                 %Coloring existence plot
stability_color = [k;zeros(steps,3)];                               %Coloring stability plot
r_stars = zeros(1,steps);

current_bifurcation_jumps = 0;
bifurcation_detected = 0;

normAs = zeros(1, steps);

disp(['Start',', norm of inverse = ' , num2str(norm((DF_arc(solutions(:,1), m, poles,u_tilde, tangents(:,1)))^-1,inf)), ', omega = ', num2str(solutions(end,1)), ', mu = ', num2str(momentum(1))])
%Continuation
i = 1;
end_segment = 0;
point_verified = 0;
while i <= steps
    %Stepping
    r_stars(i) = r_star;
    
    approx_solution = solutions(:,i) + step_size * tangents(:,i);    
    [fixed_solution, iter, normAs(i)] = newton_F_arc(approx_solution, m, poles, u_tilde, approx_solution, tangents(:,i));
    
    %Detecting bifurcations
    if current_bifurcation_jumps < max_bifurcation_jumps && detect_bifurcation == 1 && normAs(i) > 10^5 && i*step_size > 4
        [u,s,v] = svd(Df(fixed_solution, m, poles, u_tilde));
        bifurcation_detected = 1;
        current_bifurcation_jumps = current_bifurcation_jumps +1;
    end
    
    vectorized_solutions(:,1 + i*n*m: (i+1)*n*m) = vectorize(fixed_solution, m);
    
    solutions(:,i+1) = fixed_solution;
    nullspace = null(Df(fixed_solution, m, poles, u_tilde));
    
    max_dot = 0;
    max_j = 1;
    for j = 1:size(nullspace,2)
        if norm(dot(nullspace(:,j), tangents(:,i))) > max_dot
            max_dot = norm(dot(nullspace(:,j), tangents(:,i)));
            max_j = j;
        end
    end        
    tangents(:,i+1) = nullspace(:,max_j);
    
    %Dealing going along the bifurcation branch
    if bifurcation_detected == 1
        phi_1 = v(:,end);
        phi_2 = v(:,end-1);
        

        if dot(phi_2, tangents(:,i+1)) > 1e-2
            new_tangent = phi_1 - dot(phi_1,tangents(:,i+1))/dot(phi_2,tangents(:,i+1)) *phi_2;
        else
            new_tangent = phi_2 - dot(phi_2,tangents(:,i+1))/dot(phi_1,tangents(:,i+1)) *phi_1;
        end

        
        tangents(:,i+1) = new_tangent/norm(new_tangent);
        
        %load_solution = fixed_solution; load_tangent = new_tangent;
        %save p0n10m1_bif_1 load_solution load_tangent;
        
        bifurcation_detected = 0;
        %break;
    end
    
    
    if dot(tangents(:,i+1), tangents(:,i)) < 0
        tangents(:,i+1) = -tangents(:,i+1);
    end    
    
    [momentum(i+1), energy(i+1)] = plot_values(solutions(:,i+1), m, poles);  
    
    %Verifying solution and stability
    if existence_proof
        proof_data_uniform(i) = interval_pseudo_arclength2(solutions(:,i), solutions(:,i+1), m, poles, u_tilde, tangents(:,i), tangents(:,i+1), r_star);
        proof_data_point(i+1) = interval_f_U(solutions(:,i+1), m, poles, u_tilde, 1e-10);
        
        
        if 0 < proof_data_point(i+1) && proof_data_point(i+1) <= 1e-10
            [val,tmp_err] = stability_verif(solutions(:,i+1), solutions(:,i+1), m, poles, proof_data_point(i+1), 'point');
            if tmp_err == 0
                d = 10;
            end
        end
        
        if proof_data_uniform(i) > 0 && proof_data_uniform(i) <= r_star
            if stability_proof
                if point_verified
                    [values(:,i+1), err] = stability_verif(solutions(:,i), solutions(:,i+1), m, poles, proof_data_uniform(i), 'prop');
                else
                    [values(:,i+1), err] = stability_verif(solutions(:,i), solutions(:,i+1), m, poles, proof_data_uniform(i), 'segment');
                end
            else
                [values(:,i+1), err] = stability(solutions(:,i+1), m, poles);
            end
        else
            values(:,i+1) = stability(solutions(:,i+1), m, poles);
            err = 3; %numerical uncertainty
        end
    else
        [values(:,i+1), err] = stability(solutions(:,i+1), m, poles);
    end
    
    %Coloring    
    if err == 0
        k = [0 1 0];    %green   - Lyapunov stable
    elseif err == 1
        k = [1 0 0];    %red     - unstable
    elseif err == 2
        k = [1 0 1];    %magenta  - inconclusive test
    elseif err == 3
        k = [1 1 0];    %yellow - numercial uncertainty
    elseif err == 4
        k = [0 0 1];    %blue - used for testing
    end
    
    if ~existence_proof || (proof_data_uniform(i) > 0 && proof_data_uniform(i) <= r_star)
        g = k;
        stability_color(i+1,:) = k;
    else
        g = [0 0 0];
        stability_color(i+1,:) = [0 0 0];
    end
    h = repmat(g, n*m, 1);

    continuation_color(1+n*m*(i):n*m*(i+1),:) = h;
    existence_color(i+1,:) = g;
    
    
    %Adaptive
    if adaptive == 1
        skip = 0;
        
        if existence_proof && range_existence(solutions(end,i))
            if proof_data_uniform(i) <= 0
                r_star = max(3 * r_star/4, 1e-10);
                skip = 1;
            elseif proof_data_uniform(i) > r_star
                r_star = min(1.1 * r_star, 1e-5);
                skip = 1;
            else
                r_star = proof_data_uniform(i) + step_size/1000;
            end
        end


        if stability_proof && range_stability(solutions(end,i))
            if ~range_stability(solutions(end,max(1,i-1))) && ~point_verified
                if ~(0 < proof_data_point(i+1) && proof_data_point(i+1) <= 1e-10)
                    error("Trying to validate stability too close to branch switch")
                end
                
                [val,tmp_err] = stability_verif(solutions(:,i+1), solutions(:,i+1), m, poles, proof_data_point(i+1), 'point');
                
                if tmp_err == 0
                    point_verified = 1;
                else
                    error("Trying to validate stability too close to branch switch")
                end
            end
            if err == 3
                skip = 1;
            end
        else
            point_verified = 0;
        end

        if ~range_existence(solutions(end,i)) && ~range_stability(solutions(end,i))
            step_size = starting_step_size;
            if range_stability(solutions(end,max(1,i-1)))
                r_star = 3e-6;
            end
        end

        if skip
            step_size = 3 * step_size/4;
            continue
        else
            step_size = min(max_step_size, 1.1 * step_size);
            if step_size < max_step_size
                r_star = 1.1 * r_star;
            end
        end
    end
    
    
    for j = 1:length(values(:,i+1))
       if values(j,i+1) < 0 && abs(values(j,i+1)) > 1e-9
           end_segment = 1;
       end
    end
    
    if (solutions(end, i+1) < solutions(end, i))
        end_segment = 1;
    end
    
    if end_segment
        break;
    end
    
    
    disp(['Step # = ', num2str(i),', norm of inverse = ' , num2str(norm((DF_arc(solutions(:,i+1), m, poles,u_tilde, tangents(:,i+1)))^-1,inf)), ', omega = ', num2str(solutions(end,i+1)), ', mu = ', num2str(momentum(i+1))])
    if adaptive == 1
        disp(step_size)
    end
    
    i = i+1;
end

stop = i;
steps_taken{segment} = stop-1;

%Plot of the continuation on the sphere
figs(1) = figure(1);
subplot(1,number_of_segments+1,segment)
[X, Y, Z] = sphere;
surf(X,Y,Z,'FaceColor', [0 0 0], 'EdgeColor', 0.8*[1,1,1], 'FaceAlpha', 0.2);
hold on
scatter3(vectorized_solutions(1,1:n*m*stop),vectorized_solutions(2,1:n*m*stop),vectorized_solutions(3,1:n*m*stop),10, continuation_color(1:n*m*stop,:),'filled')
%scatter3(vectorized_solutions(1,1:n*m),vectorized_solutions(2,1:n*m),vectorized_solutions(3,1:n*m),50, segment_color,'filled')
scatter3(vectorized_solutions(1,1:n*m),vectorized_solutions(2,1:n*m),vectorized_solutions(3,1:n*m),50, continuation_color(1:n*m,:),'filled')
if poles == 1
    scatter3(0,0,1, 50, segment_color,'filled')
elseif poles == -1
    scatter3(0,0,1, 50, segment_color, 'filled')
elseif poles == 2
    scatter3(0,0,1, 50, segment_color,'filled')
    scatter3(0,0,-1, 50, segment_color,'filled')
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
if segment == 1
    title('Initial Segment')
elseif segment == 2
    title('First Bifurcation')
elseif segment == 3
    title('Second Bifurcation')
elseif segment == 4
    title('Third Bifurcation')
end

%{
%Plot of the error bounds along the continuation
if existence_proof
    subplot(2,max(number_of_segments, 3),max(number_of_segments, 3)+1);
    %scatter(momentum(2:end),proof_data,10,color1,'filled', 'DisplayName', 'r_star')
    scatter(momentum(2:end),proof_data,10,existence_color(2:end,:),'filled')
    hold on
    scatter(momentum(2:end),r_stars, 2, [0.5, 0.5, 0.5])    
    set(gca,'FontSize',15)
    axis square
    %axis ([0 inf 0 inf])
    set(gca,'yscale','log')
    xlabel('$$\mu$$', 'Interpreter', 'latex', 'FontSize', 25)
    ylabel('$$r_0$$', 'Interpreter', 'latex', 'FontSize', 25)
    title('Existence Bound Along the Continuation')
end



%Plot of the eigenvalues along the continuation
subplot(2,max(number_of_segments, 3),max(number_of_segments, 3)+2)
for i = 1:size(values, 1)
    %scatter(momentum, values(i,:), 5, stability_color, 'filled');
    scatter(solutions(end,:), values(i,:), 5, stability_color, 'filled');
    hold on   
end
%plot(momentum, zeros(1,steps+1), 'k')
plot(solutions(end,:), zeros(1,steps+1), 'k')
set(gca,'FontSize',15)
xlabel('$$\omega$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$\lambda$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
title('Eigenvalue-Momentum Diagram')
%}


%Plot of the energy along the continuation
subplot(1,number_of_segments+1,number_of_segments+1)
scatter(momentum(1:stop), energy(1:stop), 10, stability_color(1:stop,:), 'filled');
hold on
scatter(momentum(1), energy(1), 50, continuation_color(1,:), 'filled')
axis square
set(gca,'FontSize',15)
xlabel('$$\mu$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$H$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
title('Energy-Momentum Diagram')


if detect_bifurcation == 1
    figs(2) = figure(2);
    hold on
    scatter(momentum(2:end), normAs, 5, 'b', 'filled');
    set(gca,'yscale','log')
end

%Plot of the eigenvalues along the continuation
figs(3) = figure(3);
for i = 1:size(values, 1)
    %scatter(momentum, values(i,:), 5, stability_color, 'filled');
    scatter(solutions(end,:), values(i,:), 5, stability_color, 'filled');
    hold on   
end
%plot(momentum, zeros(1,steps+1), 'k')
plot(solutions(end,:), zeros(1,steps+1), 'k')
set(gca,'FontSize',15)
xlabel('$$\omega$$', 'Interpreter', 'latex', 'FontSize', 25)
ylabel('$$\lambda$$', 'Interpreter', 'latex', 'FontSize', 25)
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
title('Eigenvalue-Momentum Diagram')

figs(4) = figure(4);
scatter(solutions(end,:),momentum,'filled')
hold on

if save_plot && segment == number_of_segments
    savefig(figs, append(save_location, name), 'compact')
end

