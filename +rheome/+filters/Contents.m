% FILTERS  Spectral filter kernels g(lambda[,omega/t]) and their application.
%
% Static (spatial) kernels shape the eigenmode coefficients of a seeded delta to make
% an atom; joint (time-vertex) kernels shape the eigenmode x time/frequency spectrum to
% describe dynamics. Mirrors toolbox/eigen/eigfilter.
%
% Spatial kernels g(lambda):
%   rheome.filters.heat             - heat / diffusion low-pass   g(l) = exp(-t*l)
%   rheome.filters.mexhat           - Mexican-hat band-pass        g(l) = (t*l).*exp(-t*l)
%   rheome.filters.diffgauss        - difference-of-Gaussians band-pass  g(l) = exp(-t1*l) - exp(-t2*l)  (t1<t2)
%   filters.bandpass_spatial - Gaussian band-pass on wavenumber k=sqrt(l)/kmax (centre kc, width w)
% Spatial application:
%   rheome.filters.apply      - g(L) X : project -> scale -> reconstruct (single filter/bank)
%   rheome.filters.localize   - localize a filter (bank) at cortical vertices -> atoms
%
% Wavelet FRAME / filterbank (tile the spectrum, use filters as wavelets):
%   rheome.filters.frame          - design a wavelet frame  (itersine tight | mexhat | heat)
%                            COMPATIBILITY SHIM over graphfilterbank -- same output struct,
%                            same warnings. New code: rheome.graphfilterbank(lrange,'Wavelet',...)
%   rheome.filters.frame_gains    - evaluate member gains on the spectrum -> H [K x M]
%   rheome.filters.frame_bounds   - frame bounds A,B and tightness B/A (tight -> exact reconstruction)
%   rheome.filters.frame_analysis - field -> multi-scale wavelet coefficients W [rows x nT x M]
%   rheome.filters.frame_synthesis- wavelet coefficients -> reconstructed field ('tight' | 'dual')
%   rheome.filters.frame_scalogram- per-scale energy over time  E [M x nT]  (scales x time)
%
% Dirac quaternion field helpers (source-vector orientation freedom):
%   rheome.filters.tovec      - full quaternion [4nV x ..] -> physical 3-vector [3nV x ..]
%   rheome.filters.toquat     - physical 3-vector [3nV x ..] -> pure quaternion [4nV x ..] (w=0)
%   rheome.filters.steer      - right-quaternion steer: rigidly re-aim every dipole by q0
%
% Joint kernels, eigenmode x TIME ('ts')  g(lambda,t):
%   rheome.filters.diffusion  - physical-time heat (spreads)
%   rheome.filters.wave       - wave propagation (expanding wavefront)
%   rheome.filters.dampedwave - propagating + decaying wave
%   rheome.filters.kleingordon- massive (dispersive) wave with a frequency floor
%
% Joint kernels, eigenmode x FREQUENCY ('js')  g(lambda,omega):
%   rheome.filters.travwave   - traveling-wave dispersion ridge
%   rheome.filters.resonator  - damped harmonic resonator at f0 (complex)
%   rheome.filters.gabor      - spatiotemporal Gabor packet (scale x frequency)
%   rheome.filters.stmatern   - spatiotemporal 1/f (Whittle-Matern) background
%   rheome.filters.bandpass   - Gaussian temporal band-pass factor g(omega)
%
% Joint transform + application:
%   rheome.filters.jspectrum  - joint eigenmode-frequency spectrum  fft(c)/sqrt(NFFT)
%   rheome.filters.ijspectrum - inverse joint transform (-> mode-coeff time series)
%   rheome.filters.jtv        - evaluate a ts/js kernel onto the common (lambda,omega) grid
%   rheome.filters.impulse    - space-time impulse response of a joint filter at a vertex delta
%
% Analytic ground truth (unit sphere):
%   rheome.filters.sphere_kernel - closed-form Legendre-series response of a filter on the sphere
%
% Author: Diellor Basha, 2026
