function [Jx, Jy, Jz] = dif_part(Jx, Jy, Jz, vn, mode, who)
% DIF_PART  The part of an ambient field the divergence / curl is taken of: tangent, 3-D (all) or normal ("curv").
%   mode: true | "tan" (J - J_n n), false | "3d" (J), "curv" (J_n n); vn the unit area-weighted vertex normals.
    if islogical(mode) || isnumeric(mode), mode = ternary(logical(mode)); end
    mode = lower(char(mode));
    dn = vn(:,1).*Jx + vn(:,2).*Jy + vn(:,3).*Jz;                          % normal component J_n
    switch mode
        case 'tan',  Jx = Jx - dn.*vn(:,1);  Jy = Jy - dn.*vn(:,2);  Jz = Jz - dn.*vn(:,3);
        case 'curv', Jx = dn.*vn(:,1);  Jy = dn.*vn(:,2);  Jz = dn.*vn(:,3);
        case '3d'
        otherwise, error(['differential:' who ':option'], 'option must be true/false, "tan", "3d" or "curv".');
    end
end
function s = ternary(b)
    if b, s = 'tan'; else, s = '3d'; end
end
