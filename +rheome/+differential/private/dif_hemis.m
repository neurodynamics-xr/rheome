function hemis = dif_hemis(S, nV)
% DIF_HEMIS  The vertex index sets the per-hemisphere solves run over: S.Hemi, or all vertices as one.
%
%   hemis = dif_hemis(S, nV)      % {[nL x 1], [nR x 1]} or {(1:nV)'}
%
% Shared by rheome.differential.helmholtz and rheome.differential.poisson, which each pin and recentre per set.
%
% Author: Diellor Basha, 2026

    if isfield(S, 'Hemi') && ~isempty(S.Hemi)
        hemis = cellfun(@(c) double(c(:)), S.Hemi(:)', 'UniformOutput', false);
    else
        hemis = {(1:nV)'};
    end
end

% Author: Diellor Basha, 2026
