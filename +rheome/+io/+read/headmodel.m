function H = headmodel(HeadModelFile)
% IO.READ.HEADMODEL  Read a Brainstorm (unconstrained) head model / leadfield.
%
%   H = rheome.io.read.headmodel(HeadModelFile)
%
% Loads the forward model. For Dirac source mapping we need the UNCONSTRAINED
% surface leadfield: Gain [nCh x 3*nV] with three columns (x,y,z) per source vertex.
%
% OUTPUT (struct H):
%   .Gain          [nCh x 3*nV] unconstrained leadfield
%   .GridLoc       [nV x 3] source positions (cortex vertices)
%   .GridOrient    [nV x 3] source orientations if present, else []
%   .nCh, .nV      counts
%   .SurfaceFile   the cortex the grid lives on (relative Brainstorm path)
%   .HeadModelType e.g. 'surface'
%   .Comment
%
% Author: Diellor Basha, 2026

    raw = load(HeadModelFile);
    if ~isfield(raw, 'Gain') || isempty(raw.Gain)
        error('io:read:headmodel:noGain', 'No Gain matrix in %s', HeadModelFile);
    end
    G = double(raw.Gain);
    if mod(size(G, 2), 3) ~= 0
        error('io:read:headmodel:constrained', ...
            'Gain is [%d x %d] — Dirac mapping needs an UNCONSTRAINED leadfield [nCh x 3*nV].', ...
            size(G,1), size(G,2));
    end
    H = struct();
    H.Gain          = G;
    H.nCh           = size(G, 1);
    H.nV            = size(G, 2) / 3;
    H.GridLoc       = i_get(raw, 'GridLoc', []);
    H.GridOrient    = i_get(raw, 'GridOrient', []);
    H.SurfaceFile   = i_get(raw, 'SurfaceFile', '');
    H.HeadModelType = i_get(raw, 'HeadModelType', '');
    H.Comment       = i_get(raw, 'Comment', '');
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
