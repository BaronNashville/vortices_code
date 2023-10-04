poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Small omega = 0.1
%Initial Segment
poles_list{1} = 0;
n_list{1} = 10;
m_list{1} = 1;
coords_0_10_1 = [
    -0.7107   -0.6856    0.1528   -0.8321    0.9594   -0.1528    0.6002    0.5390    0.0981    0.0317
    0.3102    0.3254   -0.9880   -0.5546   -0.2769    0.9880    0.5238    0.5031   -0.4196   -0.4115
    0.6314   -0.6512    0.0229   -0.0096    0.0545    0.0215   -0.6045    0.6756   -0.9024    0.9109
    ];
u = reshape(coords_0_10_1, 30, 1);
alpha = 0;
omega = 0.1;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);

%First Bifurcation
load '../Data/p0n10m1_bif_1'
poles_list{2} = 0;
n_list{2} = 10;
m_list{2} = 1;
initial_solutions_list{2} = load_solution;
tangents_list{2} = load_tangent;

if adaptive == 0
    steps_list = {6000, 1000};
    %steps_list = {5229, 845};
else
    steps_list = {100000, 500};
    %steps_list = {30000, 5000, 5000};
end

number_of_segments = length(n_list)-1;

range_existence_list{1} = @(x) (1 <= x && x <= 2.27);
range_stability_list{1} = @(x) (1.56 <= x && x <= 2.2);
