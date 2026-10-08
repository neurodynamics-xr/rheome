function n = elements(d, rank)
% DOMAIN.ELEMENTS  How many elements of a rank a domain has.
%
%   n = rheome.domain.elements(d, 'vertex' | 'edge' | 'face' | 'element')
%
% The one arithmetic a shape check needs: a field's scalar count is this times the
% components per element of its type (rheome.fieldtype.components).
%
% Author: Diellor Basha, 2026

    % ⭐ A PRODUCT DEFERS TO ITS FACTORS. 'element' is the product count; any other rank is factor
    % 1's, so an operator built for a mesh sees exactly what it saw before and cannot tell that it is
    % now on a product. That is what makes the kind additive rather than a migration.
    if strcmp(d.kind, 'product')
        if strcmp(char(rank), 'element')
            n = d.factors{1}.nV * d.nFrames;
        else
            n = rheome.domain.elements(d.factors{1}, rank);
        end
        return
    end
    K = rheome.domain.kinds();
    r = K.ranks{K.kind == string(d.kind)};
    if ~ismember(char(rank), r)
        error('domain:elements:rank', 'A %s domain has no %s rank (it has %s).', ...
              d.kind, char(rank), strjoin(r, ', '));
    end
    switch char(rank)
        case {'vertex','element'}, n = d.nV;
        case 'edge',               n = d.nE;
        case 'face',               n = d.nF;
    end
end
% Author: Diellor Basha, 2026
