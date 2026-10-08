function S = onstate(y, opts)
% DETECT.ONSTATE  Binary on/off state of an amplitude series from the bimodality of its log.
%
%   S = rheome.detect.onstate(y)                       % y [1 x n] amplitude samples (e.g. per-cycle mean envelope)
%   S = rheome.detect.onstate(y, Posterior=0.5)
%
% An oscillation that switches on and off gives a log-amplitude distribution with two modes: the
% background floor and the rhythm. A two-component Gaussian mixture on log(y) finds them, and a sample
% is ON when its posterior for the upper component exceeds Posterior. The threshold is where the
% posterior crosses, in amplitude units.
%
% ⭐ WHY NOT A ROBUST MEDIAN THRESHOLD (rheome.detect.spindle's). median + k*MAD assumes events are SPARSE:
% the median is then background. If alpha is on for most of a recording the median sits INSIDE alpha,
% the threshold lands above the rhythm's own typical level, and a sustained episode is cut into
% sub-second fragments. The mixture places the threshold between the two states wherever they are.
%
% ⚠ NOT EVERY SERIES IS BIMODAL, and then there is no state to find. .bimodal compares the two-
% component fit against one by BIC and requires the components to be separated (Ashman's D > 2);
% when it is false the on/off labels are returned but should not be read as states.
%
% OUTPUT (struct S)
%   .on [1 x n] logical   .threshold (amplitude)   .posterior [1 x n]
%   .bimodal   .ashmanD   .dBIC (BIC one - BIC two; > 0 favours two)   .fracOn
%   .mu [lo hi] .sigma [lo hi] (of log y)   .weight [lo hi]
%
% See also: rheome.detect.occupancy, rheome.detect.spindle, alpha_occupancy_omega
%
% Author: Diellor Basha, 2026

    arguments
        y (1,:) double
        opts.Posterior (1,1) double {mustBeInRange(opts.Posterior, 0, 1)} = 0.5
        opts.Replicates (1,1) double {mustBePositive} = 3
    end
    x = log(max(y, realmin))';
    o = statset('MaxIter', 500);
    g1 = fitgmdist(x, 1);
    g2 = fitgmdist(x, 2, 'Replicates', opts.Replicates, 'Options', o, 'RegularizationValue', 1e-6);
    [mu, ord] = sort(g2.mu(:)');  sg = sqrt(squeeze(g2.Sigma(:))');  sg = sg(ord);  w = g2.ComponentProportion(ord);
    P = posterior(g2, x);  pOn = P(:, ord(2))';
    S.on = pOn > opts.Posterior;
    S.posterior = pOn;
    xs = linspace(mu(1), mu(2), 400)';  ps = posterior(g2, xs);  ps = ps(:, ord(2));
    k = find(ps > opts.Posterior, 1);  S.threshold = NaN;  if ~isempty(k), S.threshold = exp(xs(k)); end
    S.ashmanD = sqrt(2) * abs(mu(2) - mu(1)) / sqrt(sg(1)^2 + sg(2)^2);
    S.dBIC = g1.BIC - g2.BIC;
    S.bimodal = S.dBIC > 10 && S.ashmanD > 2;
    S.fracOn = mean(S.on);
    S.mu = mu;  S.sigma = sg;  S.weight = w;
end

% Author: Diellor Basha, 2026
