function [T0, sbands, sframe] = spatial(readfcn, G, g, cfg, opts)
% INGEST.SPATIAL  The wavelength axis: a tight log-itersine frame on the sensor graph, per tile.
%
%   [T0, sbands, sframe] = rheome.ingest.spatial(readfcn, G, g, cfg, Mask=[], Chunk=60)
%
% The field on the sensors at every sample is filtered by each member of a
% rheome.graphfilterbank('logitersine') on the graph's spectrum -- octaves of wavelength, tight
% (sum of squared gains = 1), with a "longer than the array" member (the low-pass, which
% holds the spatial mean) and a "shorter than the alias floor" member (the high-pass).
% Per level-0 tile and channel: the member's coefficient energy over the tile (sum) and
% the largest coefficient magnitude (max). Both merge across time and across the sensor
% tree like every other statistic.
%
% ⭐ PARSEVAL ON THE GRAPH: summed over all sensors, the members' energies equal the
% field's energy exactly (A = B = 1), so at the ROOT of the sensor tree the spatial bands
% partition sumX2 per tile. Per channel or per group they do not: a spectral filter
% moves energy between sensors, and tightness is a whole-array property. What a group's
% spatial energies say is how much of the field's energy at each wavelength its sensors
% carry.
%
% ⚠ Wavelengths in millimetres come from rheome.sensors.calibrate, which on the CTF helmet
% agrees between its two routes at only 0.36 and is accurate above ~122 mm. sbands
% carries both the dimensionless k range (exact) and the calibrated wavelengths.
%
% INPUTS:
%   readfcn  @(a, b) -> [b-a+1 x C] samples a..b, ALL the store's channels
%   G        rheome.sensors.graph on exactly those channels (in the store's channel order)
%   g        rheome.ingest.grid;  cfg  rheome.ingest.config (SpaceVoices)
%   Mask     [nT x 1] logical valid samples;  Chunk  seconds per read (60)
% OUTPUT:
%   T0      .senergy [K0 x C x nS] double, .senvMax [K0 x C x nS] single
%   sbands  table: sband_id, kind ('longer' | 'octave' | 'shorter'), k_lo, k_hi
%           (dimensionless sqrt(lambda)), k_center, wavelength_lo, wavelength_hi,
%           wavelength_center (m, NaN without calibration), n_modes
%   sframe  .A .B .alpha .agree .lambdaUsable .Q .Lambda
%
% See also: rheome.graphfilterbank, rheome.sensors.modes, rheome.sensors.calibrate, rheome.ingest.reduce
%
% Author: Diellor Basha, 2026

    arguments
        readfcn (1,1) function_handle
        G (1,1) struct
        g (1,1) struct
        cfg (1,1) struct
        opts.Mask logical = logical([])
        opts.Chunk (1,1) double {mustBePositive} = 60
    end
    nT = g.nT;  fs = g.fs;  F = g.F;  K0 = g.K0;
    C = G.nV;
    mask = opts.Mask;  if isempty(mask), mask = true(nT, 1); end
    M = rheome.sensors.modes(G);
    gfb = rheome.graphfilterbank(M.Lambda, 'Wavelet', 'logitersine', 'VoicesPerOctave', cfg.SpaceVoices);
    H = graphfilters(gfb, 'Lambda', M.Lambda);                 % [C x nS] gains on the spectrum
    Phi = M.Phi;  nS = size(H, 2);
    b = framebounds(gfb, M.Lambda);
    try
        c = rheome.sensors.calibrate(G);  alpha = c.alpha;  agree = c.agree;  lu = c.lambdaUsable;
    catch
        alpha = NaN;  agree = NaN;  lu = NaN;
    end
    % band table from the members' supports on the k axis
    kTop = sqrt(max(M.Lambda));  V = cfg.SpaceVoices;
    J = nS - 2;
    kLo = zeros(nS, 1);  kHi = zeros(nS, 1);  kind = cell(nS, 1);
    kLo(1) = 0;  kHi(1) = kTop * 2^(-(J + 0.5)/V);  kind{1} = 'longer';           % the low-pass owns below the lowest centre - half a voice
    for jj = 1:J
        m = J + 2 - jj;                                                          % member order: coarse first
        kLo(m) = kTop * 2^(-(jj + 0.5)/V);  kHi(m) = kTop * 2^(-(jj - 0.5)/V);  kind{m} = 'octave';
    end
    kLo(nS) = kTop * 2^(-0.5/V);  kHi(nS) = kTop;  kind{nS} = 'shorter';
    kC = sqrt(kLo .* kHi);  kC(1) = kHi(1) / 2;
    wl = @(k) 2*pi * sqrt(alpha) ./ max(k, eps);
    sbands = table((1:nS)', kind, kLo, kHi, kC, wl(kHi), wl(kLo), wl(kC), sum(H.^2, 1)', ...
                   'VariableNames', {'sband_id','kind','k_lo','k_hi','k_center','wavelength_lo','wavelength_hi','wavelength_center','n_modes'});
    sbands.wavelength_hi(1) = Inf;  sbands.wavelength_lo(nS) = 0;
    sframe = struct('A', b.A, 'B', b.B, 'alpha', alpha, 'agree', agree, 'lambdaUsable', lu, ...
                    'Q', qfactor(gfb), 'Lambda', M.Lambda(:)', 'voices', V, 'family', 'logitersine');
    % the transform, in time chunks aligned to level-0 frames
    T0 = struct('senergy', zeros(K0, C, nS), 'senvMax', zeros(K0, C, nS, 'single'));
    n = round(opts.Chunk * fs);  n = max(F, floor(n / F) * F);
    for a = 1:n:nT
        b2 = min(a + n - 1, nT);
        X = double(readfcn(a, b2));                                              % [len x C]
        X(~mask(a:b2), :) = 0;
        len = b2 - a + 1;  nfr = ceil(len / F);
        if nfr * F > len, X(len+1:nfr*F, :) = 0; end
        Cf = X * Phi;                                                            % [len x C] spectral coefficients
        k0 = (a - 1) / F;                                                        % first frame index - 1
        for m = 1:nS
            Y = (Cf .* H(:, m)') * Phi';                                         % member m's field, [len x C]
            E = reshape(Y.^2, F, nfr, C);
            T0.senergy(k0 + (1:nfr), :, m) = T0.senergy(k0 + (1:nfr), :, m) + reshape(sum(E, 1), nfr, C);
            T0.senvMax(k0 + (1:nfr), :, m) = max(T0.senvMax(k0 + (1:nfr), :, m), ing_bound(reshape(max(abs(reshape(Y, F, nfr, C)), [], 1), nfr, C), 'up'));
        end
    end
end
% Author: Diellor Basha, 2026
