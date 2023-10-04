function [values, err] = stability(x, m, poles)
%Function that runs the stability computations based on m

if m == 1
    [values, err] = stability_m1(x, m, poles);
elseif m == 2
    [values, err] = stability_m2(x, m, poles);
elseif m >= 3
    [values, err] = stability_m3(x, m, poles);
end