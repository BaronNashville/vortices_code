poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Triaugmented triangular prism
%Initial Segment
poles_list{1} = 0;
n_list{1} = 3;
m_list{1} = 3;
coords_0_0_3 = [
    1.0000    0.3555    0.3555   -0.5000   -0.7111   -0.7111   -0.5000    0.3555    0.3555
         0    0.6158    0.6158    0.8660    0.0000    0.0000   -0.8660   -0.6158   -0.6158
         0    0.7031   -0.7031         0    0.7031   -0.7031         0    0.7031   -0.7031
    ];
u = [coords_0_0_3(:,1);coords_0_0_3(:,2);coords_0_0_3(:,3)];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);

%First Bifurcation
load '../Data/p0n3m3_bif_1'
poles_list{2} = 0;
n_list{2} = 3;
m_list{2} = 3;
initial_solutions_list{2} = load_solution;
tangents_list{2} = load_tangent;

%Second Bifurcation
load '../Data/p0n3m3_bif_2'
poles_list{3} = 0;
n_list{3} = 3;
m_list{3} = 3;
initial_solutions_list{3} = load_solution;
tangents_list{3} = load_tangent;

if adaptive == 0
    steps_list = {9000, 1000, 10000};
else
    steps_list = {50000, 50000, 10000};
end

number_of_segments = length(n_list);

range_existence_list{1} = @(x) (0 <= x && x <= 4.98);
range_existence_list{2} = @(x) (5.03 <= x && x <= 5.17);
range_existence_list{3} = @(x) (5.225 <= x);

range_stability_list{1} = @(x) (0.3 <= x && x <= 4.98);
range_stability_list{2} = @(x) (5.03 <= x && x <= 5.17);
range_stability_list{3} = @(x) (5.225 <= x);
