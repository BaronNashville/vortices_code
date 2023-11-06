poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};


%Tetrahedron with vertex as pole
%Initial Segment
poles_list{1} = 1;
n_list{1} = 1;
m_list{1} = 3;    
coords_1_1_3 = [
    0.9428   -0.4714   -0.0000   -0.4714
         0    0.8165         0   -0.8165
   -0.3333   -0.3333    1.0000   -0.3333
   ];
u = [coords_1_1_3(:,1)];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2, 1);


if adaptive == 0
    steps_list = {5000};
else
    steps_list = {5000};
end

number_of_segments = length(n_list);

range_existence_list{1} = @(x) (0 <= x);
range_stability_list{1} = @(x) (0.1 <= x);

