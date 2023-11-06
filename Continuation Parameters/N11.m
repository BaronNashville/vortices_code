poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Triaugmented triangular prism
%Initial Segment
poles_list{1} = 1;
n_list{1} = 5;
m_list{1} = 2;
coords_1_5_2 = [
    0.4146   -0.4146   -0.9849    0.9849    0.5142   -0.5142    0.4020   -0.4020    0.5188   -0.5188    0
    0.7484   -0.7484   -0.0096    0.0096   -0.8401    0.8401    0.7257   -0.7257   -0.2874    0.2874    0
    0.5176    0.5176    0.1729    0.1729    0.1729    0.1729   -0.5583    0.5583   -0.8051   -0.8051    1
    ];
u = [coords_1_5_2(:,1);coords_1_5_2(:,3);coords_1_5_2(:,5);coords_1_5_2(:,7);coords_1_5_2(:,9)];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);

%First Bifurcation & Switching to m = 1 symmetry 
load '../Data/p1n5m2_bif_1'
poles_list{2} = 0;
n_list{2} = 11;
m_list{2} = 1;
initial_solutions_list{2} = load_solution;
tangents_list{2} = load_tangent;

if adaptive == 0
    steps_list = {3000, 1000};
else
    steps_list = {50000, 50000};
end

number_of_segments = length(n_list);

range_existence_list{1} = @(x) (0 <= x);
range_existence_list{2} = @(x) (1.32 <= x && x <= 1.47);

range_stability_list{1} = @(x) (0.5 <= x && x <= 1.26);
range_stability_list{2} = @(x) (1.35 <= x && x <= 1.46);
