function E = fold(solve, X)
% the analysis A X of fields X [V x m] on the solve's glued vertices: exact 2 P' M x; coexact G^-1 P' M x (G = P' M P on the
% coexact modes). As nxr.register.modes.fold.
P = double(solve.P); co = solve.kind(:) == 2;
E = P' * (solve.mass(:) .* X);
E(~co, :) = 2 * E(~co, :);
if any(co)
    G = P(:, co)' * (solve.mass(:) .* P(:, co)); G = (G + G') / 2;
    Rc = chol(G); E(co, :) = Rc \ (Rc' \ E(co, :));
end
end
