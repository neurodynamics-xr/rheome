function hemis = dif_hemis(S, nV)
% DIF_HEMIS  The vertex index sets the per-hemisphere solves run over: S.Hemi, or the mesh's connected components.
%
%   hemis = dif_hemis(S, nV)      % {[nL x 1], [nR x 1]}: one set per hemisphere
%
% S.Hemi is either a cell of vertex index sets (rheome.io.read.surface) or a per-vertex hemisphere index
% [nV x 1] (e.g. 1 = left, 2 = right). Without S.Hemi the sets are the mesh's connected components: the
% Laplacian's nullspace has one constant per component, so each must be pinned on its own (one solve over
% two disconnected hemispheres is singular).
% Shared by rheome.differential.helmholtz and rheome.differential.poisson, which each pin and recentre per set.
%
% Author: Diellor Basha, 2026

    if isfield(S, 'Hemi') && ~isempty(S.Hemi)
        if iscell(S.Hemi)
            hemis = cellfun(@(c) double(c(:)), S.Hemi(:)', 'UniformOutput', false);
            return
        end
        if numel(S.Hemi) ~= nV
            error('differential:hemis:size', 'S.Hemi: a cell of vertex index sets, or one hemisphere index per vertex (%d), got %d.', nV, numel(S.Hemi));
        end
        ci = S.Hemi(:);
    else
        F  = double(S.Faces);
        G  = graph([F(:,1); F(:,2); F(:,3)], [F(:,2); F(:,3); F(:,1)], [], nV);
        ci = conncomp(simplify(G))';
    end
    hemis = arrayfun(@(k) find(ci == k), unique(ci)', 'UniformOutput', false);
end

% Author: Diellor Basha, 2026
