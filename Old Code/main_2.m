%Master script that will run everything
clear, clc

existence_proof = 0;             %Set to 1 if you want to prove the existence of the equilibria
stability_proof = 0;            %Set to 1 if you want to prove the stability of the equilibria
adaptive = 0;                  %Set to 1 if you want to take adaptive step sizes to ensure proof of the above

%Note that the stability paramter cna only be set to 1 if the existence paramter is as well.
%Adaptive requires at least existence to be 1 and supports stability either way

%Setting the appropriate parameters for the continuation
%We run the file 'Continuation Parameters/N$' where $= 4,...,12 is the desired number of vertices
run('Continuation Parameters/N10')

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
    
    continuation   
end