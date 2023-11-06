This code was submitted as a part of the paper titled Determination of stable branches of relative
equilibria of the N-vortex problem on the sphere.

This READ_ME describes how to reproduce the proofs presented in the paper above.

There are 2 types of numerical proofs this code can produce.

Recall the following colour coding when looking at the plots.

    Green      - Lyapunov stable
    Red        - Unstable
    Magenta    - Inconclusive test
    Yellow     - Numerical uncertainty

If the existence of the solution could not be proven, we color the point black.

Type 1: Proving the existence of stable branches of relative equilibria paramaterized by omega for N = 4, 5, ..., 12
    To produce these proofs, we will be using the "main_continuation.m" file.

    Step 1: Select N, the number of vortices.
    Step 2: Select what level of proof you want to do, ie either only existence or both existence and stability.
    Step 3: If you want the computer to take the time needed to complete the proof everywhere, set adaptive to 1.
            Otherwise some step sizes may be too large for some segments of the branch.
    Step 4: Select whether you want to save the produced plot or not.
    Step 5: Run the file and wait for the plots to appear, you will be updated the the step sizes and value of omega
            along the way.

Type 2: Proving the existence of stable relatively equilbria for large omega for N = 10, 11, 12
    To produce these proofs, we will be using the "main_large_omega.m" file.

    Step 1: Select N, the number of vortices.
    Step 2: Select what level of proof you want to do, ie either only existence or both existence and stability.
    Step 3: Select whether you want to save the produced plot or not.
    Step 4: Run the file and wait for the plot to appear.
