function [momentum,energy] = plot_values(x, m, poles)
%Computes the momentum and energy of a given configuration of vortices
%Used for plotting

%Extracting the data from the input
n = (length(x)-2)/4;
u = reshape(x(1:3*n), 3, n);

momentum = m*sum(u(3,:));

if poles == 1
    momentum = momentum +1;
elseif poles == -1
    momentum = momentum -1;
end

if momentum < 0
    momentum = -momentum;
end

%Expressing our vdata in a form easy to work with
if poles == 0
    N = n*m;
    u_no_sym = vectorize(x,m);
elseif poles == 1
    N = n*m+1;
    u_no_sym = [vectorize(x,m),[0;0;1]];
elseif poles == -1
    N = n*m;
    u_no_sym = [vectorize(x,m),[0;0;-1]];
else
    N = n*m+2;
    u_no_sym = [vectorize(x,m),[0;0;1],[0;0;-1]];
end

energy = 0;

%Computing the energy using vectorization
for i = 1:N-1
    diff = u_no_sym(:,1) - u_no_sym;
    diff(:,1) = [];    u_no_sym(:,1) = [];
    
    norm_diff = sum(diff.^2,1).^(1/2);
    
    energy = energy + sum(log(norm_diff)/4);
end

energy = -1/2 * energy;