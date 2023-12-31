poles = 1;
n = 5;
m = 2;
coords_1_5_2 = [
    0.4146   -0.4146   -0.9849    0.9849    0.5142   -0.5142    0.4020   -0.4020    0.5188   -0.5188    0
    0.7484   -0.7484   -0.0096    0.0096   -0.8401    0.8401    0.7257   -0.7257   -0.2874    0.2874    0
    0.5176    0.5176    0.1729    0.1729    0.1729    0.1729   -0.5583    0.5583   -0.8051   -0.8051    1
    ];
u = [coords_1_5_2(:,1);coords_1_5_2(:,3);coords_1_5_2(:,5);coords_1_5_2(:,7);coords_1_5_2(:,9)];

old_input = [
    0.3898
    0.7037
    0.5941
   -0.9560
   -0.0091
    0.2932
    0.4993
   -0.8153
    0.2932
    0.4296
    0.7754
   -0.4628
    0.5834
   -0.3232
   -0.7451
    5.3788
    5.1870
    5.1870
    4.7049
    4.5249
   -0.0000
    0.6376
    ];

%Large omega = 30 - no symmetries
poles = 0;
n = 11;
m = 1;
coords_0_11_1 = [
    0.1768   -0.0308   -0.2959    0.2745    0.0440   -0.3488   -0.4329    0.4542   -0.1142   -0.0632    0.3363
    0.0327   -0.4559    0.3534   -0.3279    0.4319   -0.2793    0.0326   -0.0498    0.1363   -0.1683    0.2943
    0.9837    0.8895    0.8875    0.9040    0.9008    0.8946    0.9008    0.8895    0.9841    0.9837    0.8946
    ];
u = reshape(coords_0_11_1, 33, 1);