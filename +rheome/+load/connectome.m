function C = connectome(name)
% LOAD.CONNECTOME  Load the cached whole-brain connectome operators & eigenbases.
%
%   C = rheome.load.connectome(name)
%
% Returns the struct cached by rheome.import.connectome:
%   .meta        Kc, SmoothHops, EdgeThr, gamma, builtFrom
%   .W .keep     vertex connectome weights + largest-connected-component indices
%   .laplacian   struct(.A .B .N  .Phi .Lambda .Reproducible=true)   Connectome Laplacian
%   .lb          struct(.A .B .gamma .Phi .Lambda .Reproducible=false) LB-Connectome
%
% For the Dirac-Connectome vector basis, lift the LB-Connectome basis on demand:
%   dc = rheome.eigen.lift(C.lb.Phi, C.lb.Lambda, C.lb.B);   % quaternion basis (inherits lb's instability)
%
% See also: rheome.import.connectome, rheome.eigen.lift, rheome.filters.frame
%
% Author: Diellor Basha, 2026

    f = fullfile(rheome.load.root(), char(name), 'connectome.mat');
    if ~exist(f, 'file')
        error('load:connectome:missing', ...
            'No connectome cached for ''%s''. Run  rheome.import.connectome(''%s'')  first.', name, name);
    end
    C = getfield(builtin('load', f, 'conn'), 'conn');
end

% Author: Diellor Basha, 2026
