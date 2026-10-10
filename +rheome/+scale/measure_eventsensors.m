function [T, ev] = measure_eventsensors(name, S, Best, opts)
% SCALE.MEASURE_EVENTSENSORS  The sensor samples and channels that carry the tracked alpha events.
%
%   [T, ev] = rheome.scale.measure_eventsensors(name, S, Best)     % Best: 3rd output of measure_grouptrack
%   [T, ev] = rheome.scale.measure_eventsensors(name, S, Best, Hemi="L", Depth=8, JLevel=4, Fraction=0.5)
%
% For every configuration in Best (the highest-scoring real Viterbi path per hemisphere x depth x rate),
% rheome.detect.eventsensors on the same 8-16 Hz analytic sensor signal and the same minimum-norm kernel
% rheome.flow.envelopemodes builds the tracked envelope from (8th-order IIR band-pass, Hilbert at the
% sensors). Tracker frame f of window w at dyadic level jj covers 300-fps frames
% frame0 + (w-1)*WindowS*300 + (f-1)*2^jj + (1 : 2^jj) (frame0: grouptrack's clean-span offset), i.e. samples ((frame-1)*dec + 1) ... at the recording rate.
%
% Returns rheome.scale.rows (analysis "eventsensors"), band "<hemi>/d<depth>/<fps>fps", metrics:
%   n_steps, passed (1 if the path beat the surrogate threshold), carrying_channels_median (per step),
%   carrying_channels_union, footprint_top10_share (share of the footprint in the 10 strongest channels)
% ev: the figure data for the configuration (Hemi, Depth, JLevel) -- E (eventsensors output), the window's
% raw and band-passed sensor data, times, channel names, 2-D sensor positions (azimuthal equidistant
% about the array's mean direction), the track and its tiles. Draw with rheome.show.eventsensors.
%
% See also: rheome.detect.eventsensors, rheome.show.eventsensors, rheome.scale.measure_grouptrack
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S struct
        Best struct
        opts.Hemi (1,1) string = "L"
        opts.Depth (1,1) double = 8
        opts.JLevel (1,1) double = 4
        opts.Fraction (1,1) double = 0.5
        opts.Rate (1,1) double = 300
        opts.WindowS (1,1) double = 2
        opts.Band (1,2) double = [8 16]
    end
    st = rheome.load.study(name);  fs = st.rec.sfreq;  F = double(st.rec.F(S.iSel, :));
    lab = {st.chan.Channel(S.iSel).Name};  clear st
    Res = rheome.inverse.mne(S.G, S.ncm, struct('ChannelTypes', {S.chT}, 'InverseMeasure', 'amplitude', 'nVert', S.nV));
    bp = designfilt('bandpassiir', 'FilterOrder', 8, 'HalfPowerFrequency1', opts.Band(1), ...
                    'HalfPowerFrequency2', opts.Band(2), 'SampleRate', fs);
    Fb = filtfilt(bp, F.').';  Fa = hilbert(Fb.').';
    dec = round(fs / opts.Rate);  nW = opts.WindowS * opts.Rate;
    T = table();  ev = struct([]);
    for b = Best(:)'
        H = S.B.(char(b.hemi));  rows = reshape((double(H.gv(:))' - 1) * 3 + (1:3)', [], 1);
        w0 = b.frame0 + (b.window - 1) * nW;                                    % 300-fps frames before the window
        bs = 2^b.jj;  f0 = w0 + ((1:max(b.track.frames))' - 1) * bs;   % 300-fps frame before each tracker frame
        s0 = w0 * dec;  pad = round(0.5 * fs);
        s = max(1, s0 + 1 - pad) : min(size(F, 2), s0 + nW * dec + pad);       % the window, +-0.5 s
        fsm = [f0 * dec + 1, (f0 + bs) * dec] - (s(1) - 1);
        E = rheome.detect.eventsensors(b.track, b.G, S.G(:, rows), Res.ImagingKernel(rows, :), Fa(:, s), ...
                                       FrameSamples=fsm, Fraction=opts.Fraction);
        nc = cellfun(@numel, E.channels);  fp = sort(E.footprint, 'descend');
        T = [T; rheome.scale.rows("eventsensors", ["n_steps" "passed" "carrying_channels_median" ...
             "carrying_channels_union" "footprint_top10_share"], ...
             [numel(nc) b.passed median(nc) numel(unique(vertcat(E.channels{:}))) sum(fp(1:min(10, end))) / sum(fp)], ...
             ["steps" "bool" "channels" "channels" "fraction"], sprintf('%s/d%d/%gfps', b.hemi, b.depth, b.rate))]; %#ok<AGROW>
        if b.hemi == opts.Hemi && b.depth == opts.Depth && b.jj == opts.JLevel
            ev = struct('E', E, 'Xraw', F(:, s), 'Xband', Fb(:, s), 't', (s - 1) / fs, 'labels', {lab}, ...
                        'pos2', i_azimuthal(reshape(S.Loc, 3, [])'), 'loc', reshape(S.Loc, 3, [])', 'track', b.track, 'tiles', b.G, ...
                        'hemi', b.hemi, 'depth', b.depth, 'rate', b.rate, 'score', b.score, 'passed', b.passed);
        end
    end
end

% azimuthal equidistant projection about the mean direction of the array (a standard MEG topography layout)
function p = i_azimuthal(L)
    c = mean(L, 1);  u = L - c;  z = c / norm(c);            % mean direction, outward from the head origin
    [~, ~, V] = svd(u - (u*z') * z, 'econ');  e1 = V(:, 1)';  e2 = cross(z, e1);
    r = L ./ vecnorm(L, 2, 2);  th = acos(max(min(r * z', 1), -1));
    ph = atan2(r * e2', r * e1');  p = th .* [cos(ph) sin(ph)];
end

% Author: Diellor Basha, 2026
