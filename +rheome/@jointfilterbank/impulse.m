function U = impulse(obj, m, vertex)
% IMPULSE  Space-time impulse response of one member at a vertex.
%
%   U = impulse(jfb, m, vertex)      -> [nV x nT]
%
% Seeds a delta at the vertex with a FLAT temporal spectrum and propagates it through
% member m, giving the dynamics that member describes.
%
% ⚠ Needs a Transform (to reach the vertex domain) AND a SignalLength (so the frequency
% bins map back to time). A bank built from an arbitrary Frequencies slice has no such
% mapping and is refused rather than half-handled.
%
% See also: spectralResponse, rheome.graphtransform.eigen
%
% Author: Diellor Basha, 2026

    if isempty(obj.Transform)
        error('jointfilterbank:noTransform', ...
            ['An impulse response needs a Transform to reach the vertex domain. ' ...
             'Attach rheome.graphtransform.eigen(Phi, B, Lambda).']);
    end
    if isempty(obj.SignalLength)
        error('jointfilterbank:noGrid', ...
            ['This bank was built from an explicit Frequencies slice, so its bins do ' ...
             'not map back to a time axis. Rebuild with ''SignalLength'' to take an ' ...
             'impulse response.']);
    end
    T  = obj.Transform;
    N  = obj.SignalLength;
    nV = T.rows;
    if ~isscalar(vertex) || vertex < 1 || vertex > nV || mod(vertex,1) ~= 0
        error('jointfilterbank:vertex', 'vertex must be an integer in 1..%d.', nV);
    end

    x = zeros(nV, 1);  x(vertex) = 1;
    c0 = T.forward(x);                                   % [K x 1] seed coefficients
    W  = jointfilters(obj, m);                           % [K x nOmega]

    % Place the member on the full DFT grid, then invert along time. 'symmetric' fills
    % the negative half for a positive-half bank, which is what makes U real.
    G = zeros(numel(c0), N);
    switch obj.Half
        case 'positive', G(:, 1:obj.NumOmega) = W;
        case 'full',     G = W;
    end
    coeffs = real(ifft(G .* c0, N, 2, 'symmetric'));
    U = T.inverse(coeffs);
end

% Author: Diellor Basha, 2026
