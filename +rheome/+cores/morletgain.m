function g = morletgain(f, f0, fwhm)
% CORES.MORLETGAIN  Morlet gain per FFT bin: a Gaussian at f0 of full width fwhm, peak 2, zero at f <= 0.
%
%   g = rheome.cores.morletgain(f, f0, fwhm)
%
%   g(f) = 2 exp(-1/2 ((f - f0)/(f0/cycles))^2) (f > 0),   cycles = f0 2 sqrt(2 ln 2) / fwhm
%
% i.e. a Gaussian of standard deviation fwhm/(2 sqrt(2 ln 2)). The peak is 2 (the analytic signal doubles the positive half),
% and the gain is exactly 0 at DC and at negative frequencies. f in Hz is the FFT axis (fftfreq, not rfftfreq); f0, fwhm in Hz.
%
% Ports nxr-analysis src/nxr/vortex.py morlet_gain (line 156, commit 30b9622).
%
% See also: rheome.cores.bandframes
%
% Author: Diellor Basha, 2026
    arguments
        f double
        f0 (1,1) double {mustBePositive}
        fwhm (1,1) double {mustBePositive}
    end
    cycles = f0 * 2 * sqrt(2 * log(2)) / fwhm;
    g = 2 * exp(-0.5 * ((f - f0) / (f0 / cycles)).^2) .* (f > 0);
end
