%Master script that will run proofs for large omega steady states
clear, clc

N = 12;                          %Select N between in 10, 11, 12

existence_proof = 1;             %Set to 1 if you want to prove the existence of the equilibria
stability_proof = 1;             %Set to 1 if you want to prove the stability of the equilibria (Need to first prove existence)

save_plot = 0;                   %Set to 1 if you want to save your plot

if existence_proof && stability_proof
    save_location = append('../Large Omega Figures/', num2str(N), '/Stable Proof/');
    name = append('N', num2str(N), '_proof');
else
    save_location = append('../Large Omega Figures/', num2str(N), '/Stable Steady State/');
    name = append('N', num2str(N), '_branch');
end

%Setting the appropriate parameters for the continuation
%We run the file 'Continuation Parameters/N$' where $= 4,...,12 is the desired number of vertices
run(append('Large Omega Parameters/N', num2str(N)));

lambda = lambda_generator(u, omega, m, poles);
initial_solution = newton_f([u;lambda;alpha;omega], m, poles, u_tilde);
large_omega

