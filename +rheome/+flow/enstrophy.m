function ken = enstrophy(ctx, w)
% FLOW.ENSTROPHY  Bilinear Gram for total (or scout-weighted) enstrophy  Omega = sum_v w.*omega^2.
%   ken = rheome.flow.enstrophy(ctx [,w])   Omega(t) = ctx.F(:,t)' * ken.gram * ctx.F(:,t)
%
% omega = scalar vorticity (rheome.flow.curl). Enstrophy is quadratic in the data, so the region
% total is a single [C x C] form Q = Wcurl' diag(w) Wcurl. Per-vertex omega^2 maps come from
% ken.map on the fly.
%
% OUTPUT (struct ken): .gram [C x C], .w [V x 1], .domain='gram', .out, .chNames,
%   .map (@(F) per-vertex enstrophy map [V x nT]).
%
% See also: rheome.flow.curl, rheome.flow.energy, rheome.flow.helicity
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(w), w = i_vertexArea(ctx.fg); end
    w  = w(:);
    Wc = rheome.differential.curl(ctx.currentKernel, ctx.S, ctx.fg);    % [V x C] vorticity kernel
    Q  = Wc' * (w .* Wc);  Q = (Q + Q')/2;                % [C x C]
    ken = struct('gram', Q, 'w', w, 'domain', 'gram', 'out', 'total vorticity^2 (enstrophy)', ...
        'chNames', {ctx.chNames}, 'provenance', struct('op','sum w*omega^2'));
    ken.map = @(F) w .* (Wc * F).^2;
end

function a = i_vertexArea(fg)
    a = accumarray(fg.Faces(:), repmat(fg.FaceArea, 3, 1), [fg.nV 1]) / 3;
end

% Author: Diellor Basha, 2026
