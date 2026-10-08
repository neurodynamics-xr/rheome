function E = scaleSpectrum(obj, X)
% SCALESPECTRUM  Energy per vertex, reduced over scales.
%
%   E = scaleSpectrum(gfb, X)      X [rows x nT] -> E [nV x nT]
%
% The name transfers from cwtfilterbank VERBATIM: it averages over scale in both
% worlds, and only the surviving axis differs -- time there, vertex here. The result
% is a cortical map of where the energy is.
%
% See also: vertexSpectrum, wt
%
% Author: Diellor Basha, 2026

    T = gfb_transform(obj);
    M = obj.NumMembers;
    E = zeros(T.nV, size(X, 2));

    if isfield(T, 'filter') && ~isempty(T.filter)
        for m = 1:M
            E = E + T.norm(T.filter(obj.G_{m}, X)).^2;
        end
        return;
    end

    C = T.forward(X);
    H = graphfilters(obj, 'Lambda', gfb_lambda(obj, T));
    for m = 1:M
        E = E + T.norm(T.inverse(H(:, m) .* C)).^2;
    end
end

% Author: Diellor Basha, 2026
