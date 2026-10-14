function Xa = atlasband(X, fs, f)
% FILTERS.ATLASBAND  The analytic signal of the atlas's one-wavelet-per-octave band that holds f.
%
%   Xa = rheome.filters.atlasband(X, fs, 10)      % rows of X through the alpha wavelet (6.07-12.14 Hz)
%
% The band's member of @timefilterbank at ONE voice per octave, anchored so it is centred on the tower band
% (rheome.detect.atlastime: centre sqrt(2)/P, half-power edges 1/P and 2/P), applied to each row by FFT on the
% row's own length, zero-phase, analytic (positive frequencies only). Scaled so a sinusoid at the centre keeps
% its amplitude, as hilbert(filtfilt(bandpass)) does: abs(Xa) is the band's envelope.
% ⚠ The kernel is two octaves wide at its feet (fc/2 to 2fc) and periodic on the row: rows shorter than a few
% supports (99.9 % in 4.8 cycles, ~0.56 s at alpha) are edge-conditioned at both ends.
%
% See also: rheome.detect.atlastime, rheome.timefilterbank
%
% Author: Diellor Basha, 2026

    arguments
        X double
        fs (1,1) double {mustBePositive}
        f (1,1) double {mustBePositive}
    end
    A = rheome.detect.atlastime(f);  N = size(X, 2);
    tfb = rheome.timefilterbank(N, SamplingFrequency=fs, VoicesPerOctave=1, Anchor=A.fc);
    m = find(tfb.K_ == 0);
    if isempty(m), error('filters:atlasband:range', 'The %.3g Hz band does not fit %d samples at %g Hz.', A.fc, N, fs); end
    H = freqz(tfb);  g = zeros(1, N);  g(1:size(H, 2)) = sqrt(2) * H(m, :);
    Xa = ifft(fft(X, [], 2) .* g, [], 2);
end

% Author: Diellor Basha, 2026
