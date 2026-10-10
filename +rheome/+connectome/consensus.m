function [G, Gw] = consensus(A, D, hemi, nBins)
% CONNECTOME.CONSENSUS  Distance-dependent group consensus connectome (Betzel et al. 2019).
%
%   [G, Gw] = rheome.connectome.consensus(A, D, hemi)
%   [G, Gw] = rheome.connectome.consensus(A, D, hemi, nBins)
%
% A plain "keep edges present in >= x% of subjects" consensus over-represents short edges, since
% they are the most reproducible. Betzel's threshold keeps the group's EDGE-LENGTH DISTRIBUTION:
% edges are binned by distance, separately for intra- and inter-hemispheric pairs, and in each bin
% the consensus keeps as many edges as a subject has there on average, taking the most
% consistent ones first (ties broken by mean weight).
%
% INPUTS:
%   A      [n x n x S] per-subject connectomes (any weight; > 0 means the edge exists), symmetric
%   D      [n x n] inter-node distance (e.g. mean Euclidean centroid distance, mm)
%   hemi   [n x 1] hemisphere id per node (any two values, e.g. 0/1)
%   nBins  number of equal-width distance bins per edge class (default 40)
%
% OUTPUTS:
%   G   [n x n] logical consensus (symmetric, zero diagonal)
%   Gw  [n x n] G times the mean NON-ZERO weight across subjects
%
% Reference: Betzel, Griffa, Hagmann & Misic (2019) Netw Neurosci 3(2):475-496.
%
% Author: Diellor Basha, 2026

    if nargin < 4 || isempty(nBins), nBins = 40; end
    [n, ~, S] = size(A);
    assert(isequal(size(D), [n n]) && numel(hemi) == n, 'rheome:consensus:size', ...
        'D must be n x n and hemi n x 1 for n = %d nodes.', n);
    B = A > 0;
    consistency = mean(B, 3);
    meanWeight = sum(A .* B, 3) ./ max(sum(B, 3), 1);
    sameHemi = hemi(:) == hemi(:)';
    upper = triu(true(n), 1);
    Bflat = reshape(B, n*n, S);

    G = false(n);
    for edgeClass = [true false]                       % intra-, then inter-hemispheric
        idx = find(upper & (sameHemi == edgeClass));
        if isempty(idx), continue; end
        d = D(idx);
        bin = discretize(d, linspace(min(d), max(d) + eps(max(d)), nBins + 1));
        for b = 1:nBins
            inBin = idx(bin == b);
            if isempty(inBin), continue; end
            target = round(mean(sum(Bflat(inBin, :), 1)));
            [~, order] = sortrows([consistency(inBin) meanWeight(inBin)], [-1 -2]);
            keep = inBin(order(1:target));
            G(keep(consistency(keep) > 0)) = true;
        end
    end
    G = G | G';
    Gw = G .* meanWeight;
end

% Author: Diellor Basha, 2026
