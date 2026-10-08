function R = rateplan(cf, fs, samplesPerCycle, bandwidth)
% FLOW.RATEPLAN  Per-filter output rate for a constant-Q bank.
%
%   R = rheome.flow.rateplan(cf, fs, samplesPerCycle)
%   R = rheome.flow.rateplan(cf, fs, samplesPerCycle, bandwidth)
%
% A CWT is CONSTANT-Q, so bandwidth grows with centre frequency and a single flat rate for
% the whole bank is wasteful at the low end (the 1 Hz filter has a 0.214 Hz bandwidth and
% is otherwise carried at the full acquisition rate). The sensible rate is a fixed number
% of samples per CYCLE OF EACH FILTER'S OWN FREQUENCY:
%
%     r_m = min(samplesPerCycle * f_m,  fs)
%
% samplesPerCycle = 30 means 30 frames to render one full rotation -- enough to see a
% vortex turn or a wave cross, which is what pattern detection needs.
%
% ⚠ THE ACQUISITION RATE IS A HARD CEILING, and above it the target is simply unreachable.
% 30 samples/cycle needs f <= fs/30: that is 20 Hz at 600 Hz and 80 Hz at 2400 Hz. At
% 600 Hz a 60 Hz filter gets 10 samples per cycle, not 30, and no processing recovers it --
% it is an acquisition decision. `.capped` and `.actualPerCycle` report this rather than
% letting a silently short frame budget look like a choice.
%
% INPUTS:
%   cf               [1 x nF] centre frequencies (Hz)
%   fs               acquisition / working rate (Hz)
%   samplesPerCycle  target samples per cycle (30)
%   bandwidth        [1 x nF] half-power bandwidths (Hz), optional -- only for .criticalTotal
%
% OUTPUT (struct R):
%   .rate [1 x nF]   .decim [1 x nF] fs/rate    .capped [1 x nF] logical
%   .actualPerCycle [1 x nF]         .ceilingHz  fs/samplesPerCycle
%   .total  .flatTotal  .criticalTotal  .gain (flatTotal/total)
%
% See also: rheome.flow.phasebin, cwtfilterbank/powerbw
%
% Author: Diellor Basha, 2026

    if nargin < 4, bandwidth = []; end
    if ~isscalar(samplesPerCycle) || samplesPerCycle <= 0
        error('flow:rateplan:target', 'samplesPerCycle must be a positive scalar.');
    end
    cf = double(cf(:)).';

    R.ceilingHz      = fs / samplesPerCycle;
    want             = samplesPerCycle * cf;
    R.rate           = min(want, fs);
    R.capped         = want > fs;
    R.decim          = max(fs ./ R.rate, 1);
    R.actualPerCycle = R.rate ./ cf;
    R.total          = sum(R.rate);
    R.flatTotal      = numel(cf) * fs;
    R.gain           = R.flatTotal / R.total;
    if isempty(bandwidth), R.criticalTotal = NaN;
    else,                  R.criticalTotal = sum(double(bandwidth(:)).');
    end
end

% Author: Diellor Basha, 2026
