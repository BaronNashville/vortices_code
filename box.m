function box = box(small, large)
%Returns the coordinates of the corners of a box in the complex plane containing the interval [small, large]

small_bound = small/2;
large_bound = large + 1;

box = [small_bound-1i, large_bound-1i, large_bound+1i, small_bound+1i];
end