function out = to_input(solve, y, inputPositions)
% REGISTER.TO_INPUT  Values on the glued surface [V_glued x m] on the input vertices [V_input x m].
%
%   out = rheome.register.to_input(solve, y, inputPositions)
%
% A kept vertex takes its own value, a cut (medial wall) vertex the inverse-distance blend of its 3 nearest bridge vertices, so
% the glue, which carries a wall's rows onto the bridge, gives back y's bridge. As nxr.register.modes.to_input.
%
% See also: rheome.register.to_glued, rheome.register.solve
%
% Author: Diellor Basha, 2026
arguments, solve (1,1) struct, y (:,:) double, inputPositions (:,3) double, end
out = zeros(size(inputPositions, 1), size(y, 2));
m = solve.kept; out(solve.keep(1:m), :) = y(1:m, :);
cut = setdiff((1:size(inputPositions, 1))', solve.keep(1:m));
if ~isempty(cut) && size(y, 1) > m
    [nn, d] = knnsearch(solve.positions(m+1:end, :), inputPositions(cut, :), 'K', min(3, size(y, 1) - m));
    w = 1 ./ max(d, 1e-9); w = w ./ sum(w, 2);
    yb = y(m+1:end, :);
    for j = 1:size(nn, 2), out(cut, :) = out(cut, :) + w(:, j) .* yb(nn(:, j), :); end
end
end
