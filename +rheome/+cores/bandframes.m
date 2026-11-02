function z = bandframes(X, fs, f0, fwhm, fps)
% CORES.BANDFRAMES  The band's analytic signal at fps frames a second over the whole run (constant-Q Morlet band, exact for band-limited input).
%
%   z = rheome.cores.bandframes(X, fs, f0, fwhm, fps)
%
% X [C x n] sensor data sampled at fs Hz; the band is a Morlet of centre f0 and full width fwhm (Hz); z is [C x m] complex SINGLE,
% m = round(n*fps/fs) (half to even), frame j at t = (j-1)/fps. n*fps/fs should be an integer so that the time axis is exact.
% The run's FFT goes through rheome.cores.morletgain; bins 1 ... k-1 (0-based; DC and the Nyquist bin k are dropped),
% k = floor(min(n, m)/2), are copied into a spectrum of length m, the rest is 0, and z = ifft(.) * (m/n). A band-limited
% A cos(2 pi f0 t + phi) comes out as A exp(i (2 pi f0 t + phi)): amplitude A, not A/2. The frames are rounded to single, as
% in nxr; everything computed from them downstream (amplitude, smoothing) is compared to ~1e-5 of the peak.
%
% Ports nxr-analysis src/nxr/vortex.py band_at_frames (line 776, commit 30b9622).
%
% See also: rheome.cores.morletgain, rheome.cores.burststages
%
% Author: Diellor Basha, 2026
    arguments
        X (:,:) {mustBeNumeric}
        fs (1,1) double {mustBePositive}
        f0 (1,1) double {mustBePositive}
        fwhm (1,1) double {mustBePositive}
        fps (1,1) double {mustBePositive}
    end
    n = size(X, 2);
    m = evenround(n * fps / fs);
    kk = [0:ceil(n/2) - 1, -floor(n/2):-1];          % fftfreq(n, 1/fs)
    g = rheome.cores.morletgain(kk * fs / n, f0, fwhm);
    Y = fft(double(X), [], 2) .* g;
    k = floor(min(n, m) / 2);
    Yo = zeros(size(X, 1), m);
    Yo(:, 2:k) = Y(:, 2:k);                           % 0-based bins 1 ... k-1
    z = single(ifft(Yo, [], 2) * (m / n));
end
