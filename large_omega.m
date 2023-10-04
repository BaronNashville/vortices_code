
if existence_proof
    sol_bound = interval_f_U(initial_solution, m, poles, u_tilde, 1e-10);
else
    sol_bound = 0;
end

if stability_proof
    [values, err] = stability_verif(initial_solution, initial_solution, m, poles, sol_bound, 'point');
else
    [values, err] = stability(initial_solution, initial_solution, m, poles);    
end

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

stability_colour = error_code_to_colour(err);


[X, Y, Z] = sphere;
surf(X,Y,Z,'FaceColor', [0 0 0], 'EdgeColor', 0.8*[1,1,1], 'FaceAlpha', 0.2);
hold on
scatter3(vect(1,:), vect(2,:), vect(3,:), 5, stability_colour,'filled')

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