function T = bandsnr(name, noise)
% SCALE.BANDSNR  Per-octave SNR of the rest recording against the subject's own noise run.
%
%   T = rheome.scale.bandsnr(name)                 % noise from rheome.load.root()/<name>/noise.mat
%   T = rheome.scale.bandsnr(name, nrec)           % or a rec struct (.F .sfreq .ChannelName)
%
% The ten octaves of band_resolution_omega.m (0.125-0.25 ... 64-128 Hz). For each: Welch power
% of every good MEG channel in the band, the ratio rest/noise per channel, and the MEDIAN over
% channels, in dB of power. Channels are matched by name.
%
% ⚠ A DIFFERENT ROUTE TO THE SAME QUANTITY. The median band power over the tiles of a constant-Q
% store is another; this takes Welch band power over the whole record. Both are a median power
% ratio at the sensors; they are not bit-identical.
% ⚠ Octaves above the noise run's Nyquist are NaN, never extrapolated.
%
% Returns a table: band, fLo, fHi, snr_dB, nChannels.
%
% See also: rheome.scale.measure_bandresolution, noise_floor_omega
%
% Author: Diellor Basha, 2026

    st = rheome.load.study(name);
    if nargin < 2 || isempty(noise)
        N = builtin('load', fullfile(rheome.load.root(), name, 'noise.mat'), 'nrec');  noise = N.nrec;
    end
    isMEG = strcmpi(st.chan.Type, 'MEG');
    names = string(st.chan.Name(:)');
    good  = find(isMEG(:) & st.rec.ChannelFlag(:) == 1);
    [tf, jn] = ismember(names(good), string(noise.ChannelName));
    good = good(tf);  jn = jn(tf);
    keep = noise.ChannelFlag(jn) == 1;  good = good(keep);  jn = jn(keep);
    [ps, fs] = i_psd(st.rec.F(good,:), st.rec.sfreq);
    [pn, fn] = i_psd(noise.F(jn,:), noise.sfreq);
    lo = 0.125 * 2.^(0:9);  hi = 2*lo;
    T = table('Size', [10 5], 'VariableTypes', {'string','double','double','double','double'}, ...
              'VariableNames', {'band','fLo','fHi','snr_dB','nChannels'});
    for b = 1:10
        T.band(b) = sprintf('%g-%g Hz', lo(b), hi(b));  T.fLo(b) = lo(b);  T.fHi(b) = hi(b);
        T.nChannels(b) = numel(good);
        if hi(b) > min(st.rec.sfreq, noise.sfreq)/2, T.snr_dB(b) = NaN; continue; end
        a = mean(ps(fs >= lo(b) & fs < hi(b), :), 1);  z = mean(pn(fn >= lo(b) & fn < hi(b), :), 1);
        T.snr_dB(b) = 10*log10(median(a ./ z));
    end
end

function [P, f] = i_psd(X, fs)
% Welch, 64 s Hann windows so the 0.125-0.25 Hz octave holds 8 bins; 50 % overlap
    w = min(size(X,2), round(64*fs));
    [P, f] = pwelch(double(X)', hann(w), floor(w/2), w, fs);
end

% Author: Diellor Basha, 2026
