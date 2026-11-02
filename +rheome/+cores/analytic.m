function z = analytic(X, fs, f0, fwhm, dec)
% CORES.ANALYTIC  A band's analytic signal, decimated: ifft(fft(X) .* morletgain), every dec-th sample (complex double).
%
%   z = rheome.cores.analytic(X, fs, f0, fwhm, dec)
%
% X [C x n] sensor data at fs Hz; the band is a Morlet of centre f0 and full width fwhm (Hz); dec = 3 for nxr's 300 -> 100 Hz.
% z is [C x ceil(n/dec)] complex DOUBLE (nxr keeps complex128 here; rheome.cores.bandframes is the single-precision frame
% version that also truncates the spectrum, and agrees with this to single precision). A cos(2 pi f0 t) gives A exp(i 2 pi f0 t).
%
% Ports nxr-analysis src/nxr/vortex.py analytic (commit 3da1320; goldens tests/goldens/vortex-pair).
%
% See also: rheome.cores.morletgain, rheome.cores.bandframes, rheome.cores.pair
%
% Author: Diellor Basha, 2026
    arguments
        X (:,:) {mustBeNumeric}
        fs (1,1) double {mustBePositive}
        f0 (1,1) double {mustBePositive}
        fwhm (1,1) double {mustBePositive}
        dec (1,1) double {mustBePositive, mustBeInteger} = 3
    end
    n = size(X, 2);
    kk = [0:ceil(n/2) - 1, -floor(n/2):-1];          % fftfreq(n, 1/fs)
    g = rheome.cores.morletgain(kk * fs / n, f0, fwhm);
    z = ifft(fft(double(X), [], 2) .* g, [], 2);
    z = z(:, 1:dec:end);
end
