function s = smoothness(X, S, lbo)
% FLOW.SMOOTHNESS  How smooth a map or a vector field is on the surface, in millimetres.
%
%   s = rheome.flow.smoothness(X, S, lbo)
%       X    [nV x nT] scalar maps, or [3nV x nT] ambient vectors (rows x1,y1,z1,...)
%       S    surface: .Vertices (m) .Faces          lbo  .Mass .L (cotan stiffness) on S
%
% ⭐ .wavelengthMM IS THE DIRICHLET WAVELENGTH: 2*pi*sqrt(x'Mx / x'Lx), the wavelength of the single
% Laplace-Beltrami mode with the same ratio of gradient energy to energy (a Rayleigh quotient). A
% field that changes at the mesh spacing reads a few mm; a field smooth at the MEG point spread reads
% tens to hundreds of mm. For vectors the three ambient components are summed, so the turning of
% the arrows across folds counts as roughness -- which is what the eye sees in a quiver plot.
% .coherence (vectors only) is the magnitude-weighted mean cosine between the vectors at the two ends
% of every mesh edge: 1 = neighbouring arrows parallel, 0 = unrelated.
%
% OUTPUT (struct s): .wavelengthMM [1 x nT]   .coherence [1 x nT] (NaN for scalar maps)
%
% See also: rheome.flow.bandlimit, rheome.differential.helmholtzbands
%
% Author: Diellor Basha, 2026

    nV = size(S.Vertices, 1);  M = lbo.Mass;  L = lbo.L;
    if size(X, 1) == 3*nV
        num = 0;  den = 0;
        for c = 1:3, x = X(c:3:end, :);  num = num + sum(x .* (M*x), 1);  den = den + sum(x .* (L*x), 1); end
        E = unique(sort([S.Faces(:, [1 2]); S.Faces(:, [2 3]); S.Faces(:, [3 1])], 2), 'rows');
        a = (E(:,1) - 1) * 3;  b = (E(:,2) - 1) * 3;
        nT = size(X, 2);  s.coherence = zeros(1, nT);
        for t = 1:nT
            u = [X(a+1,t) X(a+2,t) X(a+3,t)];  v = [X(b+1,t) X(b+2,t) X(b+3,t)];
            s.coherence(t) = sum(sum(u .* v, 2)) / max(sum(vecnorm(u, 2, 2) .* vecnorm(v, 2, 2)), realmin);
        end
    elseif size(X, 1) == nV
        num = sum(X .* (M*X), 1);  den = sum(X .* (L*X), 1);  s.coherence = NaN(1, size(X, 2));
    else
        error('flow:smoothness:size', 'X has %d rows; expects nV = %d or 3nV.', size(X, 1), nV);
    end
    s.wavelengthMM = 1e3 * 2*pi * sqrt(num ./ max(den, realmin));
end

% Author: Diellor Basha, 2026
