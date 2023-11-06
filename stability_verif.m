function [values, err] = stability_verif(x_0, x_1, m, poles, bound, verif_type)
%Function that runs the rigorous stability computations based on m
%Verif types are: point, segment, or prop
%Point: Verify stability only at point x_0
%Segment: Verify along the linear segment going from x_0 to x_1
%Prop: Only verify invertibility of M_cal to carry over stabiility results

if m == 1
    [values, err] = stability_verif_m1(x_0, x_1, m, poles, bound, verif_type);
elseif m == 2
    [values, err] = stability_verif_m2(x_0, x_1, m, poles, bound, verif_type);
elseif m >= 3
    [values, err] = stability_verif_m3(x_0, x_1, m, poles, bound, verif_type);
end