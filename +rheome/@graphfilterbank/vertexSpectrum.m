function E = vertexSpectrum(obj, X)
% VERTEXSPECTRUM  Energy per scale, reduced over vertices. The timeSpectrum analogue.
%
%   E = vertexSpectrum(gfb, X)      X [rows x nT] -> E [M x nT]
%
% cwtfilterbank averages over TIME because time is its other axis; a graph filterbank's
% other axis is VERTICES. Trailing dimensions pass through, so a [rows x nT] field
% yields the [M x nT] scalogram without this class ever knowing time exists.
%
% One member is materialised at a time, so the [rows x nT x M] volume is never formed --
% safe for long records. The reduction uses Transform.norm, which is where a quaternion
% (Dirac) field becomes physical current magnitude rather than full quaternion norm.
%
% See also: scaleSpectrum, wt
%
% Author: Diellor Basha, 2026

    T = gfb_transform(obj);
    M = obj.NumMembers;
    E = zeros(M, size(X, 2));

    if isfield(T, 'filter') && ~isempty(T.filter)
        for m = 1:M
            E(m, :) = sum(T.norm(T.filter(obj.G_{m}, X)).^2, 1);
        end
        return;
    end

    C = T.forward(X);
    H = graphfilters(obj, 'Lambda', gfb_lambda(obj, T));
    for m = 1:M
        E(m, :) = sum(T.norm(T.inverse(H(:, m) .* C)).^2, 1);
    end
end

% Author: Diellor Basha, 2026
