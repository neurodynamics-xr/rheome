function [r, Nm] = rates(obj)
% RATES  The rate each member is evaluated at, without transforming anything.
%
%   [r, Nm] = rates(tfb)      r [M x 1] Hz, Nm [M x 1] samples over the record
%
% A member occupies nb bins of the record's one-sided grid, so its sub-band inverse FFT has
% length Nm = the smallest 5-smooth number >= Oversample * nb, and its coefficients arrive
% at Nm*fs/N Hz. `wt` uses exactly this; here it is available from the design alone, which
% is what a table of the bank's structure needs.
%
% ⭐ THE RATE IS PROPORTIONAL TO THE CENTRE FREQUENCY, not to the record's sampling rate:
% on the reference bank (4 voices, Oversample 2) every band member runs at 0.7*fc. That is the
% whole multirate claim -- an octave costs the same per tile wherever it sits, because the
% tile doubles exactly as the rate halves.
%
% See also: wt, bins, centerFrequencies
%
% Author: Diellor Basha, 2026

    nb = obj.Bins_(:, 2) - obj.Bins_(:, 1) + 1;
    Nm = arrayfun(@(k) tfb_smooth(ceil(obj.Oversample * k)), nb);
    r  = Nm * obj.SamplingFrequency / obj.SignalLength;
end
% Author: Diellor Basha, 2026
