function E = envelopemodes(name, opts)
% FLOW.ENVELOPEMODES  The analytic band envelope of the source estimate, as LB mode coefficients (cached).
%
%   E = rheome.flow.envelopemodes('sub01')                   % alpha, 8-16 Hz, 300 frames/s
%   E = rheome.flow.envelopemodes(name, Band=[8 16], Rate=300, Rebuild=false)
%
% Per hemisphere: band-pass and Hilbert AT THE SENSORS (exact, since the inverse is linear -- 270
% transforms instead of 61 000), the plain minimum-norm kernel, |J^| per vertex (rheome.flow.activation), and
% the projection onto the hemisphere's Laplace-Beltrami modes, chunk by chunk. The whole recording
% becomes [K x nFrames] per hemisphere, so every later pass -- data, surrogate, injection -- pools
% tiles with one [nTile x K] matrix instead of touching vertices again.
%
% CACHE: +data/<name>/alpha_envelope_modes.mat for the 8-16 Hz default (the file the alpha scripts
% already read), envelope_modes__<lo>-<hi>Hz__<rate>.mat otherwise. Rebuild=true recomputes.
%
% OUTPUT  E.L / E.R: struct(.C single [K x nFrames], .fs, .band)
%
% See also: rheome.flow.activation, rheome.inverse.mne, alpha_grouptrack_omega, alpha_jointstate_omega
%
% Author: Diellor Basha, 2026

    arguments
        name (1,1) string
        opts.Band (1,2) double = [8 16]
        opts.Rate (1,1) double {mustBePositive} = 300
        opts.Rebuild (1,1) logical = false
        opts.ChunkS (1,1) double {mustBePositive} = 2
    end
    D = fullfile(rheome.load.root(), char(name));
    if isequal(opts.Band, [8 16]) && opts.Rate == 300
        f = fullfile(D, 'alpha_envelope_modes.mat');
    else
        f = fullfile(D, sprintf('envelope_modes__%g-%gHz__%g.mat', opts.Band, opts.Rate));
    end
    if isfile(f) && ~opts.Rebuild
        S = load(f);  E = S.E;  return
    end
    B = rheome.load.bases(name);  Sf = rheome.load.surface(name);  st = rheome.load.study(name);
    isMEG = strcmpi(st.chan.Type, 'MEG');
    iSel = find(isMEG(:) & (st.rec.ChannelFlag(:) == 1));
    G = double(st.hm.Gain(iSel, :));  ncm = struct('NoiseCov', double(st.ncov.NoiseCov(iSel, iSel)));
    chT = st.chan.Type(iSel);  fs = st.rec.sfreq;  F = double(st.rec.F(iSel, :));  clear st
    Res = rheome.inverse.mne(G, ncm, struct('ChannelTypes', {chT}, 'InverseMeasure', 'amplitude', 'nVert', Sf.nV));
    Kall = Res.ImagingKernel;  clear Res G
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', opts.Band(1), ...
                    'HalfPowerFrequency2', opts.Band(2), 'SampleRate', fs);
    Fa = hilbert(filtfilt(bp, F.')).';  clear F
    Fa = Fa(:, 1:round(fs / opts.Rate):end);
    nFr = size(Fa, 2);  chunk = round(opts.ChunkS * opts.Rate);
    for hh = ["L" "R"]
        H = B.(char(hh));  gv = double(H.gv(:));  rows = reshape((gv' - 1) * 3 + (1:3)', [], 1);
        K = Kall(rows, :);  Phi = H.lbo.Phi;  Mm = H.lbo.Mass;
        C = zeros(size(Phi, 2), nFr, 'single');
        for c0 = 1:chunk:nFr
            ix = c0:min(c0 + chunk - 1, nFr);
            C(:, ix) = single(Phi' * (Mm * rheome.flow.activation(K * Fa(:, ix))));
        end
        E.(char(hh)) = struct('C', C, 'fs', opts.Rate, 'band', opts.Band);
    end
    save(f, 'E', '-v7.3');
end

% Author: Diellor Basha, 2026
