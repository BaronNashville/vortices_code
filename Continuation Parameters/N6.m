poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Octahedron (2 triangles)
poles_list{1} = 0;
n_list{1} = 2;
m_list{1} = 3;
coords_0_2_3 = [
         0    0.7071         0   -0.7071    0.7071   -0.7071
   -0.8165   -0.4082    0.8165    0.4082    0.4082   -0.4082
    0.5773   -0.5773   -0.5773    0.5773    0.5774   -0.5774
    ];
u = [coords_0_2_3(:,5);coords_0_2_3(:,2)];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);

if adaptive == 0
    steps_list = {2000};
else
    steps_list = {15000};
end

number_of_segments = length(n_list);

range_existence_list{1} = @(x) 0 <= x;
range_stability_list{1} = @(x) (0.2 <= x && x <= 1.4);
