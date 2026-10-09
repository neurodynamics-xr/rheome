function C = coefficients(B, K, F, fs, opts)
% SCALE.COEFFICIENTS  Graph-wavelet coefficient envelopes rolled up to cortical tiles: Prognome's MEG input.
%
%   C = rheome.scale.coefficients(B, K, F, fs)
%   C = rheome.scale.coefficients(B, K, F, fs, Depth=7, Band=[1 40], Rate=100)
%   S = rheome.scale.sensors(name);  st = rheome.load.study(name);
%   C = rheome.scale.coefficients(S.B, S.Res.ImagingKernel, st.rec.F(S.iSel,:), st.rec.sfreq);
%   save('sub_rheome_coeffs.mat', '-struct', 'C')     % what rheome.scale.run(..., Analyses="coefficients") writes
%
% THE CONTRACT (prognome docs/rheome-contract.md, reader prognome.data.sources.rheome.coefficients):
%   W       [nTile x nScale x nTime] single   envelope of the graph-wavelet coefficients, per tile
%   codes   [nTile x 1]   ladder codes: root 1 = cortex, 2 = left, 3 = right; children 2c, 2c+1
%   scales  [1 x nScale]  graph scales t of the bank (rheome.graphfilterbank.scales)
%   fs      scalar        sample rate of W, Hz
% plus .band .depth .widthsMM .hemisphere .bank (what the reader ignores, kept for provenance).
%
% Per hemisphere: band-pass and Hilbert AT THE SENSORS (exact, the inverse is linear), the
% unconstrained kernel K projected onto the hemisphere's Laplace-Beltrami modes per orientation, the
% gain g_m(lambda) of each bank member, back to vertices, the norm over the three orientations of the
% COMPLEX coefficient (the analytic envelope, not |J|, which pulses at twice the carrier -- see
% rheome.flow.activation), then the area-weighted mean over each tile (rheome.geom.tilemean's formula).
%
% ⚠ The envelope is taken per vertex BEFORE the tile mean. A band-pass coefficient averages toward zero
% over a tile larger than its scale (rheome.geom.tilemean); its envelope does not.
% ⚠ One bank for both hemispheres, on the smaller of their two lambda_max, so a scale means the same
% thing left and right. Members losing more than 5% of their mass beyond it (gfb.Usable) are dropped,
% and so is the lowpass (scaling) member: it has no scale t (NaN), and `scales` is a coordinate.
% ⚠ The tile ladder needs every node down to Depth to split (rheome.geom.tiles .heap): one early leaf
% would shift every code, so it is an error rather than a different answer.
%
% INPUTS
%   B    rheome.load.bases struct (.L/.R with .gv .S .lbo.Phi .lbo.Lambda .lbo.Mass .lbo.L)
%   K    [3nV x nCh] unconstrained imaging kernel, interleaved [x1;y1;z1;x2;...]
%   F    [nCh x nT] sensor data on the kernel's channels       fs  its sample rate, Hz
%   Depth  tree depth per hemisphere (7: 128 tiles each, 256 in all)
%   Band   [lo hi] Hz band-pass at the sensors before the Hilbert transform
%   Rate   output sample rate, Hz (fs/Rate must be an integer)
%   Bank   a rheome.graphfilterbank to use instead of the default mexhat bank
%
% See also: rheome.scale.run, rheome.flow.envelopemodes, rheome.geom.tree, rheome.geom.tiles, rheome.graphfilterbank
%
% Author: Diellor Basha, 2026

    arguments
        B struct
        K double
        F double
        fs (1,1) double {mustBePositive}
        opts.Depth (1,1) double {mustBeInteger, mustBePositive} = 7
        opts.Band (1,2) double = [1 40]
        opts.Rate (1,1) double {mustBePositive} = 100
        opts.ChunkS (1,1) double {mustBePositive} = 2
        opts.Bank = []
    end
    dec = fs / opts.Rate;
    if abs(dec - round(dec)) > 1e-9
        error('scale:coefficients:rate', 'fs/Rate = %g/%g is not an integer.', fs, opts.Rate);
    end
    if opts.Band(2) >= opts.Rate / 2
        error('scale:coefficients:band', 'Band edge %g Hz is not below the output Nyquist %g Hz.', opts.Band(2), opts.Rate / 2);
    end
    if size(K, 2) ~= size(F, 1)
        error('scale:coefficients:channels', 'K has %d channels, F has %d.', size(K, 2), size(F, 1));
    end
    hemis = ["L" "R"];
    gfb = opts.Bank;
    if isempty(gfb)
        gfb = rheome.graphfilterbank(min(arrayfun(@(h) max(B.(h).lbo.Lambda), hemis)));
    end
    t = scales(gfb);  w = widths(gfb);                   % sigma = sqrt(2t), metres on a cortex in metres
    keep = find(gfb.Usable & isfinite(t) & t > 0);       % wavelets only: the lowpass member has no scale (NaN)

    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', opts.Band(1), ...
                    'HalfPowerFrequency2', opts.Band(2), 'SampleRate', fs);
    Fa = hilbert(filtfilt(bp, F.')).';                     % [nCh x nT] complex, analytic
    Fa = Fa(:, 1:round(dec):end);
    nT = size(Fa, 2);  chunk = max(1, round(opts.ChunkS * opts.Rate));

    W = {};  codes = {};
    for hi = 1:2
        h = hemis(hi);  H = B.(h);
        Phi = H.lbo.Phi;  Mm = H.lbo.Mass;  lam = H.lbo.Lambda(:);
        rows = reshape((double(H.gv(:))' - 1) * 3 + (1:3)', [], 1);
        Km = cell(1, 3);
        for d = 1:3, Km{d} = Phi' * (Mm * K(rows(d:3:end), :)); end    % [nMode x nCh] per orientation
        g = zeros(numel(lam), numel(keep));
        for j = 1:numel(keep), gm = gain(gfb, keep(j)); g(:, j) = gm(lam); end

        T = rheome.geom.tree(H.S, L=H.lbo.L, M=Mm, MaxDepth=opts.Depth);
        G = rheome.geom.tiles(T, H.S, opts.Depth);
        if ~G.heap || any(G.depth ~= opts.Depth)
            error('scale:coefficients:ladder', ['%s hemisphere: the tree stopped above depth %d, so its ' ...
                  'node ids are not ladder codes. Use a smaller Depth.'], h, opts.Depth);
        end
        [code, o] = sort(G.node_id + 2.^opts.Depth * hi);        % hemisphere bit: L under 2, R under 3
        a = full(sum(Mm, 2));  P = double(G.P(:, o));
        Pw = (P ./ (P' * a)')' .* a';                            % [nTile x nV] area-weighted tile mean

        Wh = zeros(numel(code), numel(keep), nT, 'single');
        for c0 = 1:chunk:nT
            ix = c0:min(c0 + chunk - 1, nT);
            Cd = cellfun(@(k) k * Fa(:, ix), Km, 'UniformOutput', false);
            for j = 1:numel(keep)
                e = 0;
                for d = 1:3, e = e + abs(Phi * (g(:, j) .* Cd{d})).^2; end
                Wh(:, j, ix) = single(Pw * sqrt(e));
            end
        end
        W{hi} = Wh;  codes{hi} = code;                            %#ok<AGROW>
    end
    C = struct('W', cat(1, W{:}), 'codes', vertcat(codes{:}), 'scales', t(keep), 'fs', opts.Rate, ...
               'band', opts.Band, 'depth', opts.Depth + 1, 'widthsMM', 1e3 * w(keep), ...
               'hemisphere', [ones(numel(codes{1}), 1); 2 * ones(numel(codes{2}), 1)], ...
               'bank', sprintf('%s, %d voices/octave, lambda_max %g', gfb.Wavelet, gfb.VoicesPerOctave, gfb.SpectralRange));
end

% Author: Diellor Basha, 2026
