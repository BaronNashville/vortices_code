%Square at the equator
poles = 0;
n = 1;
m = 4;    
coords_0_1_4 = [
     0    -1     0     1
    -1     0     1     0
     0     0     0     0
   ];
u = [coords_0_1_4(:,1)];


%Tetrahedron with side length 2sqrt(2/3) built from the cube
poles = 0;
n = 2;
m = 2;
coords_0_2_2 = [
    (2/3)^(1/2),    0,              -(2/3)^(1/2),   0;
    0,              (2/3)^(1/2),    0,              -(2/3)^(1/2);
    3^(-1/2),       -3^(-1/2),      3^(-1/2),       -3^(-1/2)
    ];
u = [coords_0_2_2(:,1);coords_0_2_2(:,2)];


%Tetrahedron with vertex as pole
poles = 1;
n = 1;
m = 3;    
coords_1_1_3 = [
    0.9428   -0.4714   -0.0000   -0.4714
         0    0.8165         0   -0.8165
   -0.3333   -0.3333    1.0000   -0.3333
   ];
u = [coords_1_1_3(:,1)];



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