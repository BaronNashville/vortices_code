poles = 0;
n = 11;
m = 1;
coords = [
    0.176786, 0.0327058, 0.983706
    -0.0308289, -0.455902, 0.889496
    -0.295862, 0.353377, 0.887463
    0.27452, -0.327868, 0.903959
    0.043972, 0.431914, 0.900842
    -0.348784, -0.279335, 0.894607
    -0.432934, 0.0326227, 0.900835
    0.454214, -0.0497866, 0.889501
    -0.114167, 0.136306, 0.984066
    -0.0632152, -0.168303, 0.983706
    0.3363, 0.29427, 0.894599
    ]';

u = reshape(coords, 3*n,1);


alpha = 0;
omega = 50;
u_tilde = u;