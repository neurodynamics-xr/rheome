function [cut, B] = bridge_rows(solve, inputPositions)
% the input vertices the glue cut, and B [n_bridge x n_cut]: each bridge vertex's inverse-distance blend of its 3 nearest of
% them (3-D, the input positions; the Editor's glue sourcesOf; nxr.register.modes.bridge_rows)
m = solve.kept; nb = size(solve.positions, 1) - m;
cut = setdiff((1:size(inputPositions, 1))', solve.keep(1:m));
if nb == 0, B = sparse(0, numel(cut)); return, end
[nn, d] = knnsearch(inputPositions(cut, :), solve.positions(m+1:end, :), 'K', 3);
w = 1 ./ max(d, 1e-9); w = w ./ sum(w, 2);
B = sparse(repmat((1:nb)', 1, 3), nn, w, nb, numel(cut));
end
