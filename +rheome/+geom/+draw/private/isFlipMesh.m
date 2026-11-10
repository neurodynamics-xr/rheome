function tf = isFlipMesh(M)
% ISFLIPMESH  True if M is what rheome.geom.draw.flipMesh returns (private; also a flipOut info.cache of before 2026-11-10).
%
% Author: Diellor Basha, 2026
    tf = isstruct(M) && isscalar(M) && all(isfield(M, {'nV', 'nF', 'next', 'tail', 'L', 'hvV', 'angleSumI', 'signpost'}));
end
