%Triangular dipyramid (triangle at the equator)
poles = 2;
n = 1;
m = 3;
coords_2_1_3 = [
     1.0000   -0.5000   -0.5000     0    0
         0     0.8660   -0.8660     0    0  
         0         0         0      1   -1
   ];
u = [coords_2_1_3(:,1)];

%Triangular dipyramid (triangle vertex at the pole)
x_rot = 0;
y_rot = -pi/2;
poles = 1;
n = 2;
m = 2;
coords_1_2_2 = [
    0.0000   -0.0000   -0.0000   -1.0000    1.0000
         0    0.8660   -0.8660         0         0
    1.0000   -0.5000   -0.5000    0.0000   -0.0000
    ];
u = [coords_1_2_2(:,2);coords_1_2_2(:,4);];



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