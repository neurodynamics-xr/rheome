% +SPECTRAL  Periodic / aperiodic decomposition of spectra and of flow fields.
%
% Neurophysiological spectra carry a broadband aperiodic (1/f-like) background beneath any
% rhythmic peaks, so a band-pass alone does not isolate a rhythm: flow measured "in the alpha
% band" may be flow of the background. These functions parameterize the background (specparam /
% FOOOF, Donoghue et al. 2020) and use the fit to build a COMPLEMENTARY PAIR OF FILTERS that
% split a field into periodic and aperiodic parts.
%
%   aperiodic  - fit the aperiodic background of one or many power spectra
%   gains      - complementary gains h_per, h_ap from a fitted background (h_per^2 + h_ap^2 = 1)
%   decompose  - fit + split a coefficient array in one call (the usual entry point)
%   carrier    - carrier frequency estimated on the PERIODIC component only
%   synthesise - the inverse of aperiodic: real series with a prescribed 1/f background
%
% THE PARTITION IS EXACT. specparam is additive in LOG power, hence multiplicative in linear
% power, so gains formed as sqrt(peaks/total) do NOT partition. Instead the aperiodic gain comes
% from the fitted background and the periodic gain is the RESIDUAL:
%
%   h_ap  = sqrt( Pap ./ P )        h_per = sqrt( max(0, 1 - Pap./P) )
%
% giving h_per^2 + h_ap^2 = 1 pointwise and |C_per|^2 + |C_ap|^2 = |C|^2. Only the background
% must be fitted; peaks are optional and are never used by the split, which removes missed peaks
% and poor peak fits from the filter entirely.
%
% SYNTHESIS IS THE OTHER DIRECTION. rheome.spectral.synthesise takes a chi and returns series that
% aperiodic fits back to that chi, which is how a realistic background is put UNDER planted
% activity: fit the real recording per spatial mode, synthesise a matching background in that
% same mode basis (so it is in the span by construction), and add the planted field on top.
%
% ⚠ Measuring chi on source-space coefficients does not measure the source's chi. Projected
% sensor noise is white, so a poorly observed mode is dominated by it and its fitted exponent
% falls toward zero: planting chi 2.40 and inverting returns 2.40 noiseless but 1.10 at SNR 3.
% See docs/2026-09-26-feature-table-design.md section 28.
%
% APPLIED PER SPATIAL MODE (not per channel), the gains depend on lambda as well as omega, so the
% resulting filter is genuinely non-separable, and the fitted exponent becomes chi(lambda) -- the
% aperiodic exponent as a function of spatial scale.
%
% See also: rheome.filters.frame, rheome.flow.curl
%
% Author: Diellor Basha, 2026
