function kh = helicity(ctx, w)
% FLOW.HELICITY  Bilinear Gram for total (or scout-weighted) helicity  H = sum_v w.*(J . grad x J).
%   kh = rheome.flow.helicity(ctx [,w])   H(t) = ctx.F(:,t)' * kh.gram * ctx.F(:,t)
%
% Helicity density J.(grad x J) is bilinear in the data, so the region total is a [C x C] form
% Q = sym(currentKernel' diag(w (x) I3) Wcurlvec). Per-vertex helicity-density maps come from kh.map.
%
% OUTPUT (struct kh): .gram [C x C] symmetric, .w [V x 1], .domain='gram', .out, .chNames,
%   .map (@(F) per-vertex helicity map [V x nT]).
%
% See also: rheome.flow.curlvec, rheome.flow.energy, rheome.flow.enstrophy
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(w), w = i_vertexArea(ctx.fg); end
    w   = w(:);
    currentKernel = ctx.currentKernel;  Wcv = rheome.flow.curlvec(ctx);          % [3V x C] each
    w3  = repelem(w, 3);
    A   = currentKernel' * (w3 .* Wcv);                 % [C x C]  (J . curlJ), not symmetric
    Q   = (A + A')/2;
    kh = struct('gram', Q, 'w', w, 'domain', 'gram', 'out', 'total helicity J.(grad x J)', ...
        'chNames', {ctx.chNames}, 'provenance', struct('op','sum w*(J.curlJ)'));
    kh.map = @(F) w .* ( (currentKernel(1:3:end,:)*F).*(Wcv(1:3:end,:)*F) ...
                       + (currentKernel(2:3:end,:)*F).*(Wcv(2:3:end,:)*F) ...
                       + (currentKernel(3:3:end,:)*F).*(Wcv(3:3:end,:)*F) );
end

function a = i_vertexArea(fg)
    a = accumarray(fg.Faces(:), repmat(fg.FaceArea, 3, 1), [fg.nV 1]) / 3;
end

% Author: Diellor Basha, 2026
