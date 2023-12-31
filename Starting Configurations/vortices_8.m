%Octagon at the equator
poles = 0;
n = 1;
m = 8;
coords_0_1_8 = [
    1.0000    0.7071    0.0000   -0.7071   -1.0000   -0.7071   -0.0000    0.7071
         0    0.7071    1.0000    0.7071    0.0000   -0.7071   -1.0000   -0.7071
         0         0         0         0         0         0         0         0
   ];
u = [coords_0_1_8(:,1)];

%Hexagon at the equator + 2 poles
poles = 2;
n = 1;
m = 6;    
coords_2_1_6 = [
     1.0000    0.5000   -0.5000   -1.0000   -0.5000    0.5000    0    0
          0    0.8660    0.8660    0.0000   -0.8660   -0.8660    0    0
          0         0         0         0         0         0    1   -1
   ];
u = [coords_2_1_6(:,1)];


%Cube with side length 2/sqrt(3) regular orientation
u = [sqrt(2/3);0;1/sqrt(3);0;sqrt(2/3);-1/sqrt(3)];
poles = 0;
n = 2;
m = 4;
coords_0_2_4 = [
    0.0000   -0.8165   -0.8165   -0.0000   -0.0000    0.8165    0.8165    0.0000
    0.8165    0.0000    0.0000   -0.8165   -0.8165   -0.0000   -0.0000    0.8165
    0.5774   -0.5774    0.5774   -0.5774    0.5774   -0.5774    0.5774   -0.5774
    ];
u = [coords_0_2_4(:, 7);coords_0_2_4(:,8)];
   
%Cube with vertex as pole
poles = 2;
n = 2;
m = 3;
coords_2_2_3 = [
         0   -0.8165   -0.8165         0         0    0.8165    0.8165         0
   -0.0000    0.4714   -0.4714    0.0000   -0.9429    0.4714   -0.4714    0.9429
    1.0000   -0.3334    0.3334   -1.0000   -0.3333   -0.3334    0.3334    0.3333
    ];
u = [coords_2_2_3(:, 3);coords_2_2_3(:, 5)];

%Square antiprism regular orientation
poles = 0;
n = 2;
m = 4;
coords_0_2_4a = [
    0.8254    0.5836    0.0000   -0.5836   -0.8254   -0.5836   -0.0000    0.5836
         0    0.5836    0.8254    0.5836    0.0000   -0.5836   -0.8254   -0.5836
    0.5646   -0.5646    0.5646   -0.5646    0.5646   -0.5646    0.5646   -0.5646
    ];
u = [coords_0_2_4a(:,1);coords_0_2_4a(:,2)];

old_input = [
    0.7869
   -0.0000
    0.6171
    0.6122
    0.6122
   -0.5005
   14.7359
   13.4031
   -0.0000
    0.2981
    ];

%Square antiprism regular orientation 2-gon symmetry
poles = 0;
n = 4;
m = 2;
coords_0_2_4a = [
    0.8254    0.5836    0.0000   -0.5836   -0.8254   -0.5836   -0.0000    0.5836
         0    0.5836    0.8254    0.5836    0.0000   -0.5836   -0.8254   -0.5836
    0.5646   -0.5646    0.5646   -0.5646    0.5646   -0.5646    0.5646   -0.5646
    ];
u = [coords_0_2_4a(:,1);coords_0_2_4a(:,2);coords_0_2_4a(:,3);coords_0_2_4a(:,4)];

%Square antiprism edge facing top
poles = 0;
n = 4;
m = 2;
coords_0_4_2a = [
    -0.5646    0.5646   -0.5646    0.5646   -0.5646    0.5646   -0.5646    0.5646
    -0.3159    0.3158    0.7626    0.7625    0.3159   -0.3158   -0.7626   -0.7625
     0.7626    0.7625    0.3159   -0.3158   -0.7626   -0.7625   -0.3159    0.3158
    ];
u = [coords_0_4_2a(:,1);coords_0_4_2a(:,3);coords_0_4_2a(:,4);coords_0_4_2a(:,5)];
    
