poles_list = {};
n_list = {};
m_list = {};
initial_solutions_list = {};
tangents_list = {};
range_existence_list = {};
range_stability_list = {};

%Triaugmented triangular prism
%Initial Segment
poles_list{1} = 2;
n_list{1} = 1;
m_list{1} = 5;
coords_2_1_5 = [
     1.0000    0.3090   -0.8090   -0.8090    0.3090     0    0
         0     0.9511    0.5878   -0.5878   -0.9511     0    0
         0         0         0         0         0      1   -1
   ];
u = [coords_2_1_5(:,1)];
alpha = 0;
omega = 0;
u_tilde = u;
lambda = lambda_generator(u, omega, m_list{1}, poles_list{1});
initial_solutions_list{1} = [u;lambda;alpha;omega];
tangents_list{1} = zeros(4*n_list{1}+2,1);

%Bifurcation 1 & Switching to zero symmetry
load '../Data/p0n7m1_bif_1'
poles_list{2} = 0;
n_list{2} = 7;
m_list{2} = 1;

initial_solutions_list{2} = [0.3034
    0.9338
    0.1898
   -0.7947
    0.5773
    0.1875
   -0.7949
   -0.5775
    0.1861
    0.3035
   -0.9342
    0.1875
    0.9818
   -0.0000
    0.1898
   -0.0015
   -0.0011
    1.0000
    0.0024
    0.0017
   -1.0000
    3.1110
    3.1097
    3.1089
    3.1097
    3.1110
    3.5851
    2.4149
   -0.0000
    0.5851];

tangents_list{2} = 2 * [-0.0249
   -0.0571
    0.3208
   -0.0241
    0.0067
   -0.1227
   -0.0607
   -0.0441
   -0.3960
   -0.0011
   -0.0249
   -0.1225
   -0.0620
   -0.0060
    0.3209
   -0.2911
   -0.2113
   -0.0007
    0.4640
    0.3368
    0.0017
    0.1879
   -0.0716
   -0.2315
   -0.0715
    0.1880
    0.0005
    0.0001
    0.0000
    0.0009];


%Bifurcation 2
load '../Data/p0n7m1_bif_2'
poles_list{3} = 0;
n_list{3} = 7;
m_list{3} = 1;
initial_solutions_list{3} = load_solution;
tangents_list{3} = load_tangent;


%Bifurcation 3
load '../Data/p0n7m1_bif_3'
poles_list{4} = 0;
n_list{4} = 7;
m_list{4} = 1;
initial_solutions_list{4} = load_solution;
tangents_list{4} = -load_tangent;

if adaptive == 0
    steps_list = {650, 2215, 720, 100};
else
    steps_list = {5000, 1000, 1000, 1000};
end

number_of_segments = length(n_list);

range_existence_list{1} = @(x) (0 <= x);
range_existence_list{2} = @(x) (0.62 <= x && x <= 1.07) || (x >= 1.1);
range_existence_list{3} = @(x) (0.62 <= x && x <= 1.07) || (x >= 1.1);
range_existence_list{4} = @(x) (0.62 <= x && x <= 1.07) || (x >= 1.1);


range_stability_list{1} = @(x) (0.2 <= x && x <= 0.57);
range_stability_list{2} = @(x) (0.62 <= x && x <= 1.07) || (x >= 1.1);
range_stability_list{3} = @(x) (0.62 <= x && x <= 1.07) || (x >= 1.1);
range_stability_list{4} = @(x) (0.62 <= x && x <= 1.07) || (x >= 1.1);








