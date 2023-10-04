%Pentagonal dipyramid (pentagon at the equator)
poles = 2;
n = 1;
m = 5;    
coords_2_1_5 = [
     1.0000    0.3090   -0.8090   -0.8090    0.3090     0    0
         0     0.9511    0.5878   -0.5878   -0.9511     0    0
         0         0         0         0         0      1   -1
   ];
u = [coords_2_1_5(:,1)];

old_input = [
    0.9917
         0
    0.1287
   15.2525
   -0.0000
    0.3925
    ];


%Pentagonal dipyramid (pentagon vertex at the pole)
x_rot = 0;
y_rot = -pi/2;
poles = 1;
n = 3;
m = 2;    
coords_1_3_2 = [
     0.0000    0.0000   -0.0000   -0.0000    0.0000   -1.0000    1.0000
         0     0.9511    0.5878   -0.5878   -0.9511         0         0
    1.0000     0.3090   -0.8090   -0.8090    0.3090    0.0000   -0.0000
   ];
u = [coords_1_3_2(:,2);coords_1_3_2(:,3);coords_1_3_2(:,7)];



Rx = [
        1,              0,          0;
        0,              cos(x_rot), -sin(x_rot);
        0,              sin(x_rot), cos(x_rot)
        ];
        
Ry = [
        cos(y_rot),      0,         sin(y_rot);
        0,              1,          0;
        -sin(y_rot),    0,          cos(y_rot)
        ];