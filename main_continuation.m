%Master script that will run proofs for stable branches
clear, clc

N = 4;

existence_proof = 1;             %Set to 1 if you want to prove the existence of the equilibria
stability_proof = 1;             %Set to 1 if you want to prove the stability of the equilibria (Need to first prove existence)
adaptive = 1;                    %Set to 1 if you want to take adaptive step sizes to ensure proof of the above (Only useful if proving something)

save_plot = 0;                   %Set to 1 if you want to save your plot

if existence_proof && stability_proof && adaptive
    save_location = append('../Parameter Figures/', num2str(N), '/Stable Proof/');
    name = append('N', num2str(N), '_proof');
else
    save_location = append('../Parameter Figures/', num2str(N), '/Stable Branch/');
    name = append('N', num2str(N), '_branch');
end

%Setting the appropriate parameters for the continuation
%We run the file 'Continuation Parameters/N$' where $= 4,...,12 is the desired number of vertices
run(append('Continuation Parameters/N', num2str(N)));

steps_taken = {};

for segment = 1:number_of_segments
    n = n_list{segment};
    m = m_list{segment};
    poles = poles_list{segment};
    
    steps = steps_list{segment};
    
    initial_solution = initial_solutions_list{segment}; 
    load_tangent = tangents_list{segment}; 
    u_tilde = initial_solution(1:3*n);
    
    range_existence = range_existence_list{segment};
    range_stability = range_stability_list{segment};

    parameter_cont
end