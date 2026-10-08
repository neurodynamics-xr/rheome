function ke = energy(ctx, w)
% FLOW.ENERGY  Bilinear Gram for total (or scout-weighted) field energy  E = sum_v w.*|J|^2.
%   ke = rheome.flow.energy(ctx [,w])   E(t) = ctx.F(:,t)' * ke.gram * ctx.F(:,t)
%
% |J|^2 is quadratic in the data, so the whole-brain (or region) energy is a single [C x C]
% quadratic form Q = currentKernel' diag(w (x) I3) currentKernel -- precomputed once, evaluated per frame as
% F' Q F. w [V x 1] is a spatial weight (default = lumped vertex area). Per-VERTEX energy maps
% are nonlinear (a stored per-vertex Q is V*C*C ~ infeasible), so ke.map returns them on the fly.
%
% OUTPUT (struct ke): .gram [C x C] symmetric PSD, .w [V x 1], .domain='gram', .out, .chNames,
%   .map (@(F) per-vertex energy map [V x nT]).
%
% See also: rheome.flow.enstrophy, rheome.flow.helicity, rheome.flow.field
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(w), w = i_vertexArea(ctx.fg); end
    w  = w(:);
    w3 = repelem(w, 3);                              % [3V x 1] weight per component
    Q  = ctx.currentKernel' * (w3 .* ctx.currentKernel);          % [C x C]
    Q  = (Q + Q')/2;
    ke = struct('gram', Q, 'w', w, 'domain', 'gram', 'out', 'total |J|^2 energy', ...
        'chNames', {ctx.chNames}, 'provenance', struct('op','sum w|J|^2'));
    currentKernel = ctx.currentKernel;
    ke.map = @(F) w .* (currentKernel(1:3:end,:)*F).^2 + w .* (currentKernel(2:3:end,:)*F).^2 + w .* (currentKernel(3:3:end,:)*F).^2;
end

function a = i_vertexArea(fg)
    a = accumarray(fg.Faces(:), repmat(fg.FaceArea, 3, 1), [fg.nV 1]) / 3;
end

% Author: Diellor Basha, 2026
