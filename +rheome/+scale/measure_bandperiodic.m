function [T, X] = measure_bandperiodic(name, S, opts)
% SCALE.MEASURE_BANDPERIODIC  Table 5 for one subject: per band, rhythm against background and instrument.
%
%   [T, X] = rheome.scale.measure_bandperiodic(name)
%   [T, X] = rheome.scale.measure_bandperiodic(name, rheome.scale.sensors(name), Halves=true)
%   [T, X] = rheome.scale.measure_bandperiodic(name, S, Noise=nrec)     % a rec struct, as rheome.scale.bandsnr
%
% The band-SNR analysis of MS1 Table 5 with the figures removed, on the clean span periodicflow uses, and
% repeated on each half of that span for split-half reliability (MS1 group plan G4, VP1 V5). Per channel:
% Welch PSD (8 s Hann, 50%; rheome.scale.cleanwelch), aperiodic fit with a knee on [1 45] Hz minus
% LineHz and 2*LineHz +-2 Hz (rheome.spectral.aperiodic); the residual max(P - Pap, 0) IS the periodic part.
% The empty room is the subject's own noise run, Welch at the same segment length, channels matched by name.
% Per band (1-2, 2-4, 4-8, 8-16, 16-32, 32-45 Hz), each a MEDIAN OVER CHANNELS of the band-mean power:
%   total_snr_dB      total / empty room                                  dB
%   periodic_snr_dB   periodic / empty room                               dB
%   osc_snr_dB        periodic / (aperiodic + empty room) ⭐ 0 dB is the null   dB
%   periodic_fraction periodic / total                                    fraction
%   SnrFixed_total, SnrFixed_osc   10^(dB/20), the inverse's regularisation at each SNR
%   floor_total_mm, floor_osc_mm   localisation floor (rheome.inverse.resolution .lambdaLoc, left
%                     hemisphere, 150 seeds) of the minimum norm regularised at that SnrFixed   mm
%   clamped_total, clamped_osc     1 where SnrFixed fell below MinSnrFixed and was clamped to it
% Spectrum rows (band ""): aperiodic_exponent, aperiodic_offset, aperiodic_knee, fit_r2 (medians over
% channels), fit_r2_q10, n_fit_below_0p8 (channels with R^2 < 0.8), iaf_hz (peak of the median periodic
% PSD in IafRange), iaf_power_ratio (that peak's periodic / aperiodic power), n_channels, span_s.
% Analysis names: "bandperiodic" (the clean span), "bandperiodic_h1" and "bandperiodic_h2" (its halves) --
% the same metrics, so an ICC pairs h1 with h2 row for row. X is the per-(half, band) table.
%
% ⚠ THE FIT RANGE IS NOT INCIDENTAL: on a reference subject 0.5-128 Hz gives
% chi 1.1 and zero 2-4 Hz periodic power on 122/270 channels; 1-45 Hz with a knee gives chi 2.40,
% R^2 0.945 and 33. Quote chi with its range.
% ⚠ A CLAMPED FLOOR IS THE CLAMP, NOT A MEASUREMENT: two bands below MinSnrFixed report the same floor.
% The clamped_* rows say which; a group table should print them as ">= value".
% ⚠ The halves share the empty room (a separate recording) and are each refitted; each half is half
% as long, so its Welch average has half the segments.
% ⚠ For scale, a reference subject (OMEGA, 600 s): 8-16 Hz total 20.7 dB, periodic fraction 0.72,
% oscillation SNR +3.7 dB, floor 222 -> 290 mm; every other band periodic fraction <= 0.39, osc < 0.
%
% See also: rheome.scale.measure_periodicflow, rheome.scale.bandsnr, rheome.spectral.aperiodic,
%           rheome.inverse.resolution, rheome.scale.cleanwelch
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Noise = []
        opts.Bands (:,2) double = [1 2; 2 4; 4 8; 8 16; 16 32; 32 45]
        opts.FitRange (1,2) double = [1 45]
        opts.LineHz (1,1) double = 60
        opts.SegS (1,1) double = 8
        opts.IafRange (1,2) double = [7 14]
        opts.MinSnrFixed (1,1) double = 0.2
        opts.NumSeeds (1,1) double = 150
        opts.Halves (1,1) logical = true
        opts.Floors (1,1) logical = true
        opts.EdgeS (1,1) double = 10
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    st = rheome.load.study(name);  fs = st.rec.sfreq;
    names = string(st.chan.Name(S.iSel));  F = double(st.rec.F(S.iSel, :));  clear st
    nz = opts.Noise;
    if isempty(nz)
        N = builtin('load', fullfile(rheome.load.root(), name, 'noise.mat'), 'nrec');  nz = N.nrec;
    end
    [tf, jn] = ismember(names, string(nz.ChannelName));
    if ~all(tf), error('scale:bandperiodic:channels', '%d of %d channels are not in the noise run', sum(~tf), numel(tf)); end
    nf = 2^nextpow2(opts.SegS * fs);
    nfe = 2^nextpow2(opts.SegS * nz.sfreq);
    Fe = double(nz.F(jn, :));
    if size(Fe, 2) < nfe, error('scale:bandperiodic:noise', 'noise run of %.1f s is shorter than one %g s segment', size(Fe,2)/nz.sfreq, opts.SegS); end
    [Pe, fe] = pwelch(Fe.', hann(nfe), nfe/2, nfe, nz.sfreq);  clear Fe

    CS = rheome.scale.cleanspan(F, fs, EdgeS=opts.EdgeS);
    spans = {CS.first:CS.last};  an = "bandperiodic";
    if opts.Halves
        mid = floor((CS.first + CS.last) / 2);
        spans = [spans {CS.first:mid, mid+1:CS.last}];  an = [an "bandperiodic_h1" "bandperiodic_h2"];
    end
    SL = S.B.L.S;  floorCache = containers.Map('KeyType', 'double', 'ValueType', 'double');
    T = rheome.scale.rows("bandperiodic", strings(0,1), [], "");  X = table();
    for h = 1:numel(spans)
        [P, f] = rheome.scale.cleanwelch(F, fs, spans{h}, CS.mask, nf);
        keep = f >= opts.FitRange(1) & f <= opts.FitRange(2);
        for l = opts.LineHz * [1 2], keep = keep & ~(f > l-2 & f < l+2); end
        ap = rheome.spectral.aperiodic(P(keep, :), f(keep), struct('knee', true));
        fk = f(keep);  Pk = P(keep, :);  Pap = ap.Pap;  Pper = max(Pk - Pap, 0);
        nb = size(opts.Bands, 1);  R = table();
        for i = 1:nb
            lo = opts.Bands(i,1);  hi = opts.Bands(i,2);
            m = fk >= lo & fk < hi;  me = fe >= lo & fe < hi;
            tot = median(mean(Pk(m,:), 1));  per = median(mean(Pper(m,:), 1));
            apw = median(mean(Pap(m,:), 1));  noi = median(mean(Pe(me,:), 1));
            % ⚠ realmin, never eps: MEG power is ~1e-27 T^2/Hz and an eps guard returns the guard
            r = struct('band', string(sprintf('%g-%g Hz', lo, hi)), 'total_snr_dB', 10*log10(tot/noi), ...
                       'periodic_snr_dB', 10*log10(max(per, realmin)/noi), ...
                       'osc_snr_dB', 10*log10(max(per, realmin)/(apw + noi)), 'periodic_fraction', per/tot);
            r.SnrFixed_total = 10^(r.total_snr_dB/20);  r.SnrFixed_osc = 10^(r.osc_snr_dB/20);
            [r.floor_total_mm, r.clamped_total] = i_floor(r.SnrFixed_total);
            [r.floor_osc_mm,   r.clamped_osc]   = i_floor(r.SnrFixed_osc);
            R = [R; struct2table(r)]; %#ok<AGROW>
        end
        mp = median(Pper, 2);  ia = fk >= opts.IafRange(1) & fk <= opts.IafRange(2);  fi = fk(ia);
        [pk, k] = max(mp(ia));  iaf = NaN;  ipr = NaN;
        if pk > 0, iaf = fi(k);  ipr = pk / median(Pap(find(ia, 1) - 1 + k, :)); end
        sm = ["aperiodic_exponent" "aperiodic_offset" "aperiodic_knee" "fit_r2" "fit_r2_q10" ...
              "n_fit_below_0p8" "iaf_hz" "iaf_power_ratio" "n_channels" "span_s"];
        sv = [median(ap.exponent) median(ap.offset) median(ap.knee) median(ap.r2) prctile(ap.r2, 10) ...
              sum(ap.r2 < 0.8) iaf ipr numel(names) numel(spans{h})/fs];
        su = ["chi" "log10 power" "Hz^chi" "R2" "R2" "channels" "Hz" "ratio" "channels" "s"];
        T = [T; rheome.scale.rows(an(h), sm, sv, su)]; %#ok<AGROW>
        bm = ["total_snr_dB" "periodic_snr_dB" "osc_snr_dB" "periodic_fraction" "SnrFixed_total" ...
              "SnrFixed_osc" "floor_total_mm" "floor_osc_mm" "clamped_total" "clamped_osc"];
        bu = ["dB" "dB" "dB" "fraction" "ratio" "ratio" "mm" "mm" "flag" "flag"];
        for i = 1:nb
            T = [T; rheome.scale.rows(an(h), bm, R{i, bm}, bu, R.band(i))]; %#ok<AGROW>
        end
        X = [X; [table(repmat(an(h), nb, 1), 'VariableNames', {'analysis'}) R]]; %#ok<AGROW>
        fprintf('[bandperiodic %s] %-16s chi %.2f R2 %.3f IAF %.2f Hz | 8-16 Hz pf %.2f osc %+.1f dB\n', name, an(h), ...
            median(ap.exponent), median(ap.r2), iaf, R.periodic_fraction(R.band == "8-16 Hz"), R.osc_snr_dB(R.band == "8-16 Hz"));
    end

    function [mm, cl] = i_floor(sf)
    % the localisation floor at one SnrFixed (memoised: the same SnrFixed gives the same inverse)
        cl = double(sf < opts.MinSnrFixed);  sf = max(sf, opts.MinSnrFixed);  mm = NaN;
        if ~opts.Floors || ~isfinite(sf), return; end
        if isKey(floorCache, sf), mm = floorCache(sf); return; end
        Res = rheome.inverse.mne(S.G, S.ncm, struct('ChannelTypes', {S.chT}, 'InverseMeasure', 'amplitude', ...
                                                    'nVert', S.nV, 'SnrFixed', sf));
        o = rheome.inverse.resolution(Res.ImagingKernel, S.G, S.B.L.lbo, Vertices=SL.Vertices, Faces=SL.Faces, ...
                Normals=SL.VertNormals, GlobalIdx=S.B.L.gv, NumSeeds=opts.NumSeeds, Modes=0);
        mm = 1e3 * o.lambdaLoc;  floorCache(sf) = mm;
    end
end

% Author: Diellor Basha, 2026
