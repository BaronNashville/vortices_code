poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Triangular dipyramid (triangle vertex at the pole)
%Initial Segment
poles_list{1} = 1;
n_list{1} = 2;
m_list{1} = 2;
coords_1_2_2 = [
    0.0000   -0.0000   -0.0000   -1.0000    1.0000
         0    0.8660   -0.8660         0         0
    1.0000   -0.5000   -0.5000    0.0000   -0.0000
    ];
u = [coords_1_2_2(:,2);coords_1_2_2(:,4);];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);

%First Bifurcation
load '../Data/p1n2m2_bif2'
poles_list{2} = 1;
n_list{2} = 1;
m_list{2} = 4;
initial_solutions_list{2} = load_solution;
tangents_list{2} = load_tangent;

if adaptive == 0
    steps_list = {2000, 2000};
else
    steps_list = {5000, 2000};
end

range_existence_list{1} = @(x) (0 <= x && x <= 0.499);
range_existence_list{2} = @(x) (0.5 <= x);

range_stability_list{1} = @(x) (0.2 <= x && x <= 0.35) || (0.355 <= x && x <= 0.47);
range_stability_list{2} = @(x) (0.51 <= x);

range_stability_list{1} = @(x) (0.1 <= x && x <= 0.49);

number_of_segments = length(n_list);