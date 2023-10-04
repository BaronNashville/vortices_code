function [values, err] = stability_verif(x_0, x_1, m, poles, bound, verif_type)
if m == 1
    [values, err] = stability_verif_m1(x_0, x_1, m, poles, bound, verif_type);
elseif m == 2
    [values, err] = stability_verif_m2(x_0, x_1, m, poles, bound, verif_type);
elseif m >= 3
    [values, err] = stability_verif_m3(x_0, x_1, m, poles, bound, verif_type);
end