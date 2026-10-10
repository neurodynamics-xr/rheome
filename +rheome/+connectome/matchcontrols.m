function [sets, nCandidates] = matchcontrols(L, N, targets, pool, nSets, tol, seed)
% CONNECTOME.MATCHCONTROLS  Control edges matched to target edges on length and streamline-count decile.
%
%   [sets, nCandidates] = rheome.connectome.matchcontrols(L, N, targets, pool, nSets)
%   [sets, nCandidates] = rheome.connectome.matchcontrols(L, N, targets, pool, nSets, tol, seed)
%
% For H3/TP3 (control tracts): each target edge (e.g. EC -> j) gets, in each of nSets sets, one
% control edge drawn at random from the POOL edges (e.g. both terminals in Braak V-VI regions) whose
% mean streamline length is within +/- tol of the target's and whose streamline count falls in the
% same decile. Deciles are over all existing edges of N (upper triangle, N > 0).
%
% INPUTS:
%   L        [n x n] mean streamline length per edge (tck2connectome -scale_length -stat_edge mean)
%   N        [n x n] streamline count per edge
%   targets  [k x 2] node pairs of the target edges
%   pool     [m x 2] node pairs the controls are drawn from
%   nSets    number of control sets (the cards ask 100)
%   tol      relative length tolerance (default 0.10)
%   seed     rng seed (default 0)
%
% OUTPUTS:
%   sets         [k x nSets] row of POOL matched to each target in each set; NaN if no candidate
%   nCandidates  [k x 1] number of eligible pool edges per target (report it: few means weak matching)
%
% Author: Diellor Basha, 2026

    if nargin < 6 || isempty(tol), tol = 0.10; end
    if nargin < 7 || isempty(seed), seed = 0; end
    stream = RandStream('mt19937ar', 'Seed', seed);
    existing = N(triu(true(size(N)), 1));
    q = quantile(existing(existing > 0), 0.1:0.1:0.9);
    decile = @(x) 1 + sum(x(:) > q, 2);
    at = @(M, P) M(sub2ind(size(M), P(:, 1), P(:, 2)));

    tLen = at(L, targets);  tDec = decile(at(N, targets));
    pLen = at(L, pool);     pDec = decile(at(N, pool));   pOk = at(N, pool) > 0;
    k = size(targets, 1);
    sets = nan(k, nSets);  nCandidates = zeros(k, 1);
    for t = 1:k
        cand = find(pOk & abs(pLen - tLen(t)) <= tol * tLen(t) & pDec == tDec(t));
        nCandidates(t) = numel(cand);
        if isempty(cand) || ~(at(N, targets(t, :)) > 0), continue; end
        sets(t, :) = cand(randi(stream, numel(cand), 1, nSets));
    end
end

% Author: Diellor Basha, 2026
