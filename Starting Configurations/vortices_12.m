%12-gon at the equator
poles = 0;
n = 1;
m = 12;
coords_0_1_12 = [
    1.0000    0.8660    0.5000    0.0000   -0.5000   -0.8660   -1.0000   -0.8660   -0.5000   -0.0000    0.5000    0.8660
         0    0.5000    0.8660    1.0000    0.8660    0.5000    0.0000   -0.5000   -0.8660   -1.0000   -0.8660   -0.5000
         0         0         0         0         0         0         0         0         0         0         0         0
    ];
u = [coords_0_1_12(:,1)];

%10-gon at the equator + 2 poles
poles = 2;
n = 1;
m = 10;
coords_2_1_10 = [
    1.0000    0.8090    0.3090   -0.3090   -0.8090   -1.0000   -0.8090   -0.3090    0.3090    0.8090    0    0
         0    0.5878    0.9511    0.9511    0.5878         0   -0.5878   -0.9511   -0.9511   -0.5878    0    0
         0         0         0         0         0         0         0         0         0         0    1   -1
    ];
u = [coords_2_1_10(:,1)];

%Icosahedron regular orientation (vertex as pole)
poles = 2;
n = 2;
m = 5;
coords_2_2_5 = [
    0.8944   -0.8944    0.2764   -0.2764   -0.7236    0.7236   -0.7236    0.7236    0.2764   -0.2764    0    0
         0         0    0.8507   -0.8507    0.5257   -0.5257   -0.5257    0.5257   -0.8507    0.8507    0    0
    0.4472   -0.4472    0.4472   -0.4472    0.4472   -0.4472    0.4472   -0.4472    0.4472   -0.4472    1   -1
    ];
u = [coords_2_2_5(:, 1);coords_2_2_5(:,2)];

%Icosahedron (edge as pole)
poles = 0;
n = 6;
m = 2;
coords_0_6_2 = [
    -0.4472    0.4472   -0.4472    0.4472   -0.4472    0.4472   -0.4472    0.4472   -0.4472    0.4472   -1.0000    1.0000
    -0.2764    0.2764    0.7236   -0.7236    0.7236   -0.7236   -0.2763    0.2763   -0.8945    0.8945         0         0
     0.8506   -0.8506    0.5258   -0.5258   -0.5257    0.5257   -0.8506    0.8506   -0.0000    0.0000    0.0000   -0.0000
    ];
u = [coords_0_6_2(:,1);coords_0_6_2(:,2);coords_0_6_2(:,3);coords_0_6_2(:,4);coords_0_6_2(:,9);coords_0_6_2(:,11)];

%Icosahedron (face as pole)
poles = 0;
n = 4;
m = 3;
coords_0_4_3 = [
    0.3035   -0.3035    0.3036   -0.3036   -0.4911    0.4911   -0.9822    0.9822   -0.4912    0.4912   -0.6071    0.6071
   -0.5258    0.5258    0.5257   -0.5257    0.8506   -0.8506    0.0001   -0.0001   -0.8507    0.8507         0         0
    0.7946   -0.7946    0.7947   -0.7947    0.1876   -0.1876   -0.1876    0.1876    0.1875   -0.1875    0.7946   -0.7946
    ];
u = [coords_0_4_3(:, 1);coords_0_4_3(:,2);coords_0_4_3(:,5);coords_0_4_3(:,6)];

%Large omega = 40 - n3m4p0 symmetry
poles = 0;
n = 4;
m = 3;
coords_0_4_3 = [
    0.0388    0.1118   -0.1505   -0.2765    0.3724   -0.0959    0.2381   -0.4101   -0.0474   -0.3306    0.1720    0.3780
    0.1514   -0.1093   -0.0422    0.2703    0.1043   -0.3746    0.3361    0.0381    0.4091   -0.2456   -0.3742   -0.1635
    0.9877    0.9877    0.9877    0.9222    0.9222    0.9222    0.9113    0.9113    0.9113    0.9113    0.9113    0.9113
    ];
u = [coords_0_4_3(:,1);coords_0_4_3(:,4);coords_0_4_3(:,7);coords_0_4_3(:,10)];
