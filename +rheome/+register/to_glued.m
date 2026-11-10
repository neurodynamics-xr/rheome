function y = to_glued(solve, x, inputPositions)
% REGISTER.TO_GLUED  Values on the input vertices [V_input x m] on the solve's glued surface [V_glued x m].
%
%   y = rheome.register.to_glued(solve, x, inputPositions)
%
% A kept vertex takes its own value, a bridge vertex the blend of its 3 nearest cut walls (the glue's sourcesOf).
% As nxr.register.modes.to_glued.
%
% See also: rheome.register.to_input, rheome.register.solve
%
% Author: Diellor Basha, 2026
arguments, solve (1,1) struct, x (:,:) double, inputPositions (:,3) double, end
[cut, B] = bridge_rows(solve, double(single(inputPositions)));
y = [x(solve.keep(1:solve.kept), :); B * x(cut, :)];
end
