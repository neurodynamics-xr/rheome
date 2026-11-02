function [frame, sense, centre, area] = patches(C, thr, A, vertArea, Vmm)
% CORES.PATCHES  Cores of a curl movie: connected same-sense patches above threshold (nxr _core_patches).
%
%   [frame, sense, centre, area] = rheome.cores.patches(C, thr, A, vertArea, Vmm)
%
%   C [V x T] curl (single), thr threshold, A [V x V] adjacency (rheome.cores.adjacency), vertArea [V x 1] (mm^2),
%   Vmm [V x 3] positions (mm).
%   frame [n x 1] (column of C, 1-based), sense [n x 1] (+1 counter-clockwise, -1 clockwise), centre [n x 3] mm,
%   area [n x 1] mm^2.
%
% A vertex joins its above-threshold neighbours in the SAME frame; one component search per sense, +1 first. Weight =
% area*|curl|; centre = sum(w*Vmm)/sum(w) (a 3-D mean); area = sum of vertex areas. Output order: by frame, +1 before -1,
% and within a frame and sense by component label order (first vertex appearance): this sets greedy-linking priority.
%
% Ports nxr-analysis src/nxr/vortex.py _core_patches (commit 903fdc2).
%
% See also: rheome.cores.stats, rheome.cores.threshold
%
% Author: Diellor Basha, 2026
    C = single(C); thr = single(thr); vertArea = double(vertArea(:));
    fr = cell(2,1); sn = cell(2,1); ce = cell(2,1); ar = cell(2,1);
    [ei, ej] = find(triu(A, 1));
    for q = 1:2
        sg = 3 - 2*q;                                   % +1 then -1
        on = (sg * C) >= thr;                           % [V x T]
        [v, f] = find(on);                              % column-major = frame, then vertex
        if isempty(v), continue, end
        node = zeros(size(C)); node(sub2ind(size(C), v, f)) = 1:numel(v);
        both = on(ei, :) & on(ej, :);
        [e, ff] = find(both);
        G = graph(node(sub2ind(size(C), ei(e), ff)), node(sub2ind(size(C), ej(e), ff)), [], numel(v));
        lab = conncomp(G).';
        [~, ~, lab] = unique(lab, 'stable');            % labels in order of first appearance, as scipy
        nc = max(lab);
        w = vertArea(v) .* abs(double(C(sub2ind(size(C), v, f))));
        W = accumarray(lab, w, [nc 1]);
        cen = zeros(nc, 3);
        for d = 1:3, cen(:, d) = accumarray(lab, w .* Vmm(v, d), [nc 1]) ./ W; end
        fr{q} = accumarray(lab, f, [nc 1], @max);  sn{q} = repmat(sg, nc, 1);
        ce{q} = cen;  ar{q} = accumarray(lab, vertArea(v), [nc 1]);
    end
    frame = vertcat(fr{:}); sense = vertcat(sn{:}); centre = vertcat(ce{:}); area = vertcat(ar{:});
    if isempty(frame), frame = zeros(0,1); sense = zeros(0,1); centre = zeros(0,3); area = zeros(0,1); return, end
    [~, o] = sortrows([frame, -sense]);                 % stable: ties keep label order
    frame = frame(o); sense = sense(o); centre = centre(o, :); area = area(o);
end
