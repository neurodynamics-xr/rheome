function pl = sen_plane(arr)
% The array's OWN plane: two orthonormal in-array directions and a centred coordinate.
%
% ⚠ DO NOT ROTATE A LATTICE INTO ITS SVD BASIS. A square grid has DEGENERATE singular
% values, so V's in-plane columns are an arbitrary rotation and the array renders as a
% diamond -- the same degeneracy that makes an SVD projection unusable in sen_faces. For a
% 2D array, drop the coordinate axis most aligned with the plane normal and keep the other
% two raw. A chain has no degeneracy, so its first singular vector IS the shank axis.
    Pc = arr.Pos - mean(arr.Pos, 1);
    [~,~,V] = svd(Pc, 'econ');
    if arr.Dim < 2
        e1 = V(:,1);  e2 = V(:,2);
    else
        [~, drop] = max(abs(V(:,3)));
        keepAx = setdiff(1:3, drop);
        E = eye(3);  e1 = E(:,keepAx(1));  e2 = E(:,keepAx(2));
    end
    pl = struct('e1', e1, 'e2', e2, 'x1', Pc*e1, 'x2', Pc*e2, 'Pc', Pc);
end

% Author: Diellor Basha, 2026
