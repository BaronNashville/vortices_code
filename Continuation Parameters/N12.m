poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Icosahedron (face as pole)
%Initial segment
poles_list{1} = 0;
n_list{1} = 4;
m_list{1} = 3;
coords_0_4_3 = [
    0.3035   -0.3035    0.3036   -0.3036   -0.4911    0.4911   -0.9822    0.9822   -0.4912    0.4912   -0.6071    0.6071
   -0.5258    0.5258    0.5257   -0.5257    0.8506   -0.8506    0.0001   -0.0001   -0.8507    0.8507         0         0
    0.7946   -0.7946    0.7947   -0.7947    0.1876   -0.1876   -0.1876    0.1876    0.1875   -0.1875    0.7946   -0.7946
    ];
u = [coords_0_4_3(:, 1);coords_0_4_3(:,2);coords_0_4_3(:,5);coords_0_4_3(:,6)];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);

%First Bifurcation & Switching to m = 1 symmetry
load '../Data/p0n12m1_bif_1_2'
poles_list{2} = 0;
n_list{2} = 12;
m_list{2} = 1;
initial_solutions_list{2} = load_solution;
tangents_list{2} = load_tangent;


if adaptive == 0
    steps_list = {5000, 100};
else
    steps_list = {40000, 1000};
end

number_of_segments = length(n_list)-1;

range_existence_list{1} = @(x) (0 <= x);

range_stability_list{1} = @(x) (0.2 <= x && x <= 0.200465) || (0.3 <= x && x <= 0.43)|| (0.47 <= x && x <= 1.98) || (2.82 <= x);
range_stability_list{1} = @(x) (0.9 <= x && x <= 2.84);


