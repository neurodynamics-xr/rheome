function [w0, info] = carrier(Cper, f)
% SPECTRAL.CARRIER  Carrier frequency, estimated on the PERIODIC component only.
%
%   [w0, info] = rheome.spectral.carrier(Cper, f)
%
% Returns the power-weighted mean frequency of Cper. Demodulation shifts the frequency origin to
% w0, so a biased w0 biases every envelope downstream.
%
% WHY THE PERIODIC COMPONENT -- ROBUSTNESS, NOT BIAS MAGNITUDE. Measured on a synthetic 1/f^1.5
% background with a 10 Hz peak (spectral_validate, section 3):
%
%   band (Hz)    w0 periodic   w0 total    bias
%   8-13            10.183      10.185    +0.002     <- negligible at a narrow band
%   7-14            10.068      10.007    -0.061
%   5-16            10.027       9.514    -0.513
%   2-25            10.027       7.973    -2.054
%   1-40            10.027       7.064    -2.962
%
% The periodic estimate is BAND-INDEPENDENT (10.03 across every band); the total-power estimate
% slides from 10.19 to 7.06 and is therefore an artefact of where the band edges were put. The
% direction of the bias is set by whether the aperiodic centroid of the retained band sits above
% or below the peak -- it is downward for wide bands, and vanishes for a narrow one where the two
% happen to coincide. Use the periodic component because it does not depend on the band, not
% because the bias is large at any particular band.
%
% INPUTS:
%   Cper  [K x nF] periodic component (from rheome.spectral.decompose)
%   f     [nF x 1] frequencies (Hz)
%
% OUTPUT:
%   w0    scalar, Hz
%   info  .w0total  the same estimate on |Cper|^2 + |Cap|^2 if a second argument set is supplied
%         .power    [1 x nF] summed power over modes
%         .halfCycle  500/w0 (ms)
%
% See also: rheome.spectral.decompose
%
% Author: Diellor Basha, 2026

    f  = double(f(:)).';
    pw = sum(abs(Cper).^2, 1);
    if ~(sum(pw) > 0)
        error('spectral:carrier:empty', 'periodic component carries no power.');
    end
    w0 = sum(pw .* f) / sum(pw);
    info.power     = pw;
    info.halfCycle = 500 / max(w0, eps);
end

% Author: Diellor Basha, 2026
