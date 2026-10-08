function T = groups(T0, tree, chanRows)
% INGEST.GROUPS  Roll level statistics up the sensor tree: one column per internal node.
%
%   Tg = rheome.ingest.groups(T, tree, chanRows)    T = {1 x Lmax+1} from rheome.ingest.rollup (or one level)
%
% For every internal node of the tree and every level, sums and maxima over the node's
% member channels: sumX, sumX2, energy add; absMax, max, envMax take the max; min the
% min. n and nCoi are per record and are shared. chanRows maps the store's channel axis
% (1..C) to rows of the array the tree was built on (meta.iChannel against the array's
% channel list), so a tree built on the helmet's 270 sensors serves a store of any
% channel subset; nodes with no member in the store get no column.
%
% ⭐ Exact by the same argument as the time roll-up: a parent node equals its children
% merged, because every stored statistic is a sum, a count, a max or a min.
%
% OUTPUT: Tg {1 x Lmax+1} (or one struct), fields sumX, sumX2 [K x nG] double;
%         absMax, min, max [K x nG] single; energy [K x nG x B] double; envMax single;
%         .nodes [1 x nG] the internal node ids in column order
%
% See also: rheome.sensors.tree, rheome.ingest.rollup, rheome.ingest.build
%
% Author: Diellor Basha, 2026

    single_ = ~iscell(T0);
    if single_, T0 = {T0}; end
    internal = find(~tree.is_leaf)';
    nG = numel(internal);
    C = numel(chanRows);
    M = zeros(nG, C);                                   % membership: node x store channel
    for i = 1:nG
        mem = tree.members{internal(i)};
        M(i, :) = ismember(chanRows(:)', mem);
    end
    keep = any(M, 2);                                   % a node with no member in the store has no row
    internal = internal(keep);  M = M(keep, :);  nG = numel(internal);
    Ms = sparse(M);
    T = cell(size(T0));
    for L = 1:numel(T0)
        P = T0{L};  K = size(P.sumX, 1);  B = size(P.energy, 3);
        Q = struct();
        Q.sumX  = P.sumX  * Ms';
        Q.sumX2 = P.sumX2 * Ms';
        for f = ["sumAbs","sumSqrt","sumX3","sumX4"]            % sums over members, like sumX
            if isfield(P, f), Q.(f) = P.(f) * Ms'; end
        end
        Q.energy = i_sum3(P.energy, Ms);
        Q.absMax = i_max2(P.absMax, M, 0);
        Q.max    = i_max2(P.max,    M, -Inf);
        Q.min    = -i_max2(-P.min,  M, -Inf);
        Q.envMax = i_max3(P.envMax, M);
        if isfield(P, 'senergy')
            Q.senergy = i_sum3(P.senergy, Ms);
            Q.senvMax = i_max3(P.senvMax, M);
        end
        Q.n = P.n;  Q.nCoi = P.nCoi;  Q.level = P.level;  Q.K = P.K;
        Q.nodes = internal;
        T{L} = Q;
    end
    if single_, T = T{1}; end
end

function S = i_sum3(A, Ms)
    [K, C, B] = size(A);
    S = reshape(reshape(permute(A, [1 3 2]), K * B, C) * Ms', K, B, []);
    S = permute(S, [1 3 2]);
end

function R = i_max2(A, M, neutral)
    K = size(A, 1);  nG = size(M, 1);
    R = repmat(cast(neutral, 'like', A), K, nG);
    for i = 1:nG
        R(:, i) = max(A(:, M(i, :) > 0), [], 2);
    end
end

function R = i_max3(A, M)
    [K, ~, B] = size(A);  nG = size(M, 1);
    R = zeros(K, nG, B, 'like', A);
    for i = 1:nG
        R(:, i, :) = max(A(:, M(i, :) > 0, :), [], 2);
    end
end
% Author: Diellor Basha, 2026
