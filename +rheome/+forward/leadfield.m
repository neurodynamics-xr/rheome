function [G, idx, names] = leadfield(study, varargin)
% FORWARD.LEADFIELD  The MEG gain matrix for the good channels: the forward direction of source mapping.
%
%   [G, idx, names] = rheome.forward.leadfield(study)
%   [G, idx, names] = rheome.forward.leadfield(study, GlobalVertices=gv)   % one hemisphere's columns
%
% ⭐ THIS EXISTS BECAUSE THE REGISTRY DEMANDED IT. `rheome.operators.registry` requires every arrow to name a
% resolvable function, and the leadfield -- ambientVertexWorld -> sensorScalar, the forward direction
% of source mapping -- had no function to point at: it was reached as `study.hm.Gain` inline at each
% call site, and the whole forward direction was therefore missing from the arrow table while every
% inverse was present. Naming it makes the chain checkable.
%
% ⚠ The columns are [x1 y1 z1 x2 y2 z2 ...] per vertex, unconstrained, and the units are tesla per
% A m. A source vector J [3nV x nT] in A m gives B = G*J in tesla.
%
% INPUTS
%   study            rheome.load.study output (.hm.Gain, .chan.Type, .rec.ChannelFlag)
%   GlobalVertices   (optional) vertex indices to keep, e.g. rheome.load.bases hemisphere .gv
%
% OUTPUTS
%   G      [nCh x 3nV] or [nCh x 3numel(gv)]   idx  the channel rows kept   names  their names
%
% See also: rheome.forward.diracgain, rheome.forward.simulate, rheome.operators.registry, rheome.inverse.mne
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('GlobalVertices', [], @isnumeric);
    p.parse(varargin{:});
    ok  = strcmpi(study.chan.Type, 'MEG') & (study.rec.ChannelFlag(:)' == 1);
    idx = find(ok);
    G   = double(study.hm.Gain(idx, :));
    gv  = p.Results.GlobalVertices;
    if ~isempty(gv)
        gv = double(gv(:))';
        G = G(:, reshape((gv-1)*3 + (1:3)', 1, []));
    end
    names = string(study.chan.Name(idx));
end

% Author: Diellor Basha, 2026
