function W = wt(obj, X)
% WT  Graph wavelet transform: a field -> its multi-scale coefficients.
%
%   W = wt(gfb, X)      X [rows x nT] -> W [rows x nT x M]
%
% Projects once, then filters by each member:
%   C = forward(X);   W(:,:,m) = inverse(g_m(lambda) .* C)
%
% Trailing dimensions pass through untouched -- time, trials and subjects are
% properties of the SIGNAL, not axes of this transform. A graph filterbank's own
% axes are (vertex, scale).
%
% See also: iwt, vertexSpectrum, scaleSpectrum, graphfilters
%
% Author: Diellor Basha, 2026

    T = gfb_transform(obj);

    % A polynomial transform has no coefficient space, so it filters directly.
    if isfield(T, 'filter') && ~isempty(T.filter)
        M = obj.NumMembers;
        W = zeros(size(X,1), size(X,2), M);
        if ~isreal(X), W = complex(W); end
        for m = 1:M
            W(:, :, m) = T.filter(obj.G_{m}, X);
        end
        return;
    end

    C = T.forward(X);                                   % [K x nT]
    H = graphfilters(obj, 'Lambda', gfb_lambda(obj, T));% [K x M]
    M  = obj.NumMembers;  nT = size(X, 2);
    W  = zeros(size(X,1), nT, M);
    if ~isreal(X) || ~isreal(C), W = complex(W); end
    for m = 1:M
        W(:, :, m) = T.inverse(H(:, m) .* C);
    end
end

% Author: Diellor Basha, 2026
