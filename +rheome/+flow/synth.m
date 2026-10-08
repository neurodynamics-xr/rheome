function [c, tOut] = synth(C, f, fsOut, dur, opts)
% FLOW.SYNTH  Return a joint spectrum to the time axis at ANY output rate.
%
%   [c, t] = rheome.flow.synth(C, f, fsOut, dur)
%   [c, t] = rheome.flow.synth(C, f, fsOut, dur, opts)
%
% Evaluates the inverse transform of a band-selected (and optionally demodulated) joint spectrum
% on a chosen time grid. The output rate is a PARAMETER, not a preprocessing commitment: the
% band-limiting already happened in frequency, so no aliasing is introduced by a coarse grid.
%
% WHY THE RATE IS FREE. After demodulation the content sits at baseband -- a few Hz wide -- so
% Nyquist applies to the ENVELOPE bandwidth, not to the carrier. Sampling an alpha envelope at the
% acquisition rate oversamples it by two orders of magnitude, and the surplus is position jitter
% with no signal in it (the tracking floor goes as v_floor ~ edgeLength * fsOut).
%
% INPUTS:
%   C      [K x nOmega] complex joint spectrum (from rheome.flow.joint, optionally via rheome.flow.demod)
%   f      [1 x nOmega] bin frequencies (Hz) -- may be negative after demodulation
%   fsOut  output sampling rate (Hz)
%   dur    output duration (s)
%   opts   .t0  start time (default 0)
%
% OUTPUT:
%   c      [K x nT'] complex coefficients on the output grid; abs(c) is the envelope
%   tOut   [1 x nT'] output times (s)
%
% Evaluates the sum directly rather than via ifft, because the retained bins are sparse and the
% output grid is arbitrary (and generally not commensurate with the original transform length).
%
% See also: rheome.flow.joint, rheome.flow.demod
%
% Author: Diellor Basha, 2026

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'t0') || isempty(opts.t0), opts.t0 = 0; end
    f = double(f(:)).';
    if size(C,2) ~= numel(f)
        error('flow:synth:size', 'C has %d columns but f has %d entries.', size(C,2), numel(f));
    end
    nOut = max(1, round(dur * fsOut));
    tOut = opts.t0 + (0:nOut-1)/fsOut;
    c    = C * exp(1i*2*pi*(f.' * tOut));      % [K x nOmega] * [nOmega x nT']
end

% Author: Diellor Basha, 2026
