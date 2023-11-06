poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Square antiprism regular orientation
%Initial Segment
poles_list{1} = 0;
n_list{1} = 2;
m_list{1} = 4;
coords_0_2_4a = [
    0.8254    0.5836    0.0000   -0.5836   -0.8254   -0.5836   -0.0000    0.5836
         0    0.5836    0.8254    0.5836    0.0000   -0.5836   -0.8254   -0.5836
    0.5646   -0.5646    0.5646   -0.5646    0.5646   -0.5646    0.5646   -0.5646
    ];
u = [coords_0_2_4a(:,1);coords_0_2_4a(:,2)];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);


%First Bifurcation & Switching to rings of 2 vortices
load '../Data/p0n4m2_bif_1'
poles_list{2} = 0;
n_list{2} = 4;
m_list{2} = 2;
initial_solutions_list{2} = load_solution;
tangents_list{2} = load_tangent;

if adaptive == 0
    steps_list =  {2200, 720};
else
    steps_list = {10000, 20000};
end

number_of_segments = length(n_list);

range_existence_list{1} = @(x) (0 <= x);
range_existence_list{2} = @(x) (1.615 <= x);

range_stability_list{1} = @(x) (0.5 <= x && x <= 1.6);
range_stability_list{2} = @(x) (1.615 <= x && x <= 1.6151) || (1.722 <= x && x <= 1.775) || (1.818 <= x && x <= 1.9);

range_stability_list{1} = @(x) (0.2 <= x && x <= 1.6);
range_stability_list{2} = @(x) (1.62 <= x && x <= 1.9);













