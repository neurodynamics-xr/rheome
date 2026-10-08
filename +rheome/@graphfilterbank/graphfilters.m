function H = graphfilters(obj, varargin)
% GRAPHFILTERS  The bank's filters, in either domain.
%
%   H = graphfilters(gfb)                          spectral gains [K x M]
%   H = graphfilters(gfb, 'Lambda', lam)           on a supplied spectrum
%   H = graphfilters(gfb, 'Dual', true)            canonical dual  g_m / S
%   A = graphfilters(gfb, 'Type','vertex', 'Vertex', v)    the atoms [nV x P x M]
%
% 'Type' names the DOMAIN: 'spectral' (default) is the gain against lambda; 'vertex' is
% g_m(L) applied to a delta at each seed.
%
% ⚠ THERE IS NO cwtfilters. R2023b's cwtfilterbank splits the two domains across two
% accessors -- freqz for the frequency response, wavelets for the time-domain wavelet.
% This method's original "one unified accessor" rationale was written from a misremembered
% name; see docs/2026-08-22-graphfilterbank-design.md section 4.3. impulse is this class's
% counterpart to wavelets. The 'vertex' mode stays because it seeds SEVERAL vertices at
% once, which wavelets has no need to do.
%
% ⚠ 'vertex' NEEDS A LOCATION. A CWT is translation-invariant, so one wavelet per
% scale describes the filter everywhere. A graph is not: g_m(L)'s impulse response
% genuinely differs at every vertex, and there is no canonical centring.
%
% See also: impulse, gain, spectralResponse, framebounds
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Lambda', []);
    p.addParameter('Dual',   false);
    p.addParameter('Type',   'spectral');
    p.addParameter('Vertex', []);
    p.parse(varargin{:});
    o = p.Results;

    if ~any(strcmpi(o.Type, {'spectral','vertex'}))
        error('graphfilterbank:type', ...
            'Type must be ''spectral'' or ''vertex'', got ''%s''.', o.Type);
    end
    if strcmpi(o.Type, 'vertex')
        H = i_vertexfilters(obj, o);
        return;
    end

    lam = o.Lambda;
    if isempty(lam)
        if obj.HasSpectrum
            lam = obj.Lambda_;
        else
            lam = linspace(0, obj.Lmax_, 512)';   % dense grid: no spectrum was supplied
        end
    end
    lam = double(lam(:));

    M = obj.NumMembers;
    H = zeros(numel(lam), M);
    for m = 1:M
        v = obj.G_{m}(lam);
        H(:, m) = v(:);
    end

    if o.Dual
        S    = sum(H.^2, 2);
        bad  = S <= 1e-12 * max([S; realmin]);
        Sinv = zeros(size(S));
        Sinv(~bad) = 1 ./ S(~bad);
        H = H .* Sinv;
    end
end

function A = i_vertexfilters(obj, o)
% The atoms: g_m(L) applied to a delta at each seed. The cwtfilterbank/wavelets
% analogue -- but a graph has no translation invariance, so the impulse response
% genuinely differs at every vertex and a location is REQUIRED. There is no
% canonical centring either.
    T = gfb_transform(obj);
    if isempty(o.Vertex)
        error('graphfilterbank:vertex', ...
            ['''Type'',''vertex'' needs a ''Vertex''. Unlike a CWT, a graph has no ' ...
             'translation invariance, so there is no single wavelet per scale -- the ' ...
             'impulse response differs at every vertex.']);
    end
    v  = o.Vertex(:);
    nV = T.rows;
    if any(v < 1 | v > nV | mod(v,1) ~= 0)
        error('graphfilterbank:vertex', 'Vertex must be integer indices in 1..%d.', nV);
    end

    P = numel(v);
    X = full(sparse(v, (1:P)', 1, nV, P));      % a delta at each seed
    M = obj.NumMembers;

    % ⚠ THE CHEBYSHEV ROUTE HAS NO SPECTRAL DOMAIN: identity forward/inverse, all the work
    % in .filter. Multiplying a gain vector into an untransformed delta is silently wrong,
    % so dispatch on the transform rather than assuming one of the two.
    if isfield(T, 'filter') && ~isempty(T.filter)
        A = zeros(nV, P, M);
        for m = 1:M
            A(:, :, m) = T.filter(obj.G_{m}, X);
        end
        return;
    end

    lam = o.Lambda;
    if isempty(lam), lam = gfb_lambda(obj, T); end
    H = graphfilters(obj, 'Lambda', lam, 'Dual', o.Dual);
    C = T.forward(X);                            % [K x P]
    A = zeros(nV, P, M);
    if ~isreal(C), A = complex(A); end
    for m = 1:M
        A(:, :, m) = T.inverse(H(:, m) .* C);
    end
end

% Author: Diellor Basha, 2026
