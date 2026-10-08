function v = vortex(J, S, props, op)
% DETECT.VORTEX  Physical properties of the vortices in a source vector field.
%
%   v = rheome.detect.vortex(J, S)             % all properties
%   v = rheome.detect.vortex(J, S, props)      % only the requested properties
%   v = rheome.detect.vortex(J, S, props, op)  % reuse a precomputed rheome.detect.operator(S)
%
% Finds the vortex critical points (rheome.detect.criticalPoints -- the Dirac-connection route) and
% MEASURES each one's physical properties on the cortical surface, in real units:
%   position   [nc x 3]  core location (m, ambient)
%   size       [nc x 1]  characteristic radius in MM: geodesic distance from the core at which
%                        the tangential flow speed peaks (= sigma of a Gaussian-stream vortex)
%   chirality  [nc x 1]  rotation sense (r x v).n_outward: +1 counter-clockwise, -1 clockwise
%                        (viewed from OUTSIDE the surface; robust to mesh face-winding)
%   strength   [nc x 1]  circulation magnitude (|curl| at the core, field units / m)
%
% props selects which to compute/return: a char ('size'), cellstr ({'size','chirality'}), or
% 'all' (default = every property). 'position' is always returned (the anchor). Case-insensitive.
%
% The sphere (rheome.import.sphere) is the CALIBRATION STANDARD: plant an analytic vortex of known
% sigma and handedness and confirm size/chirality read back correctly (see detect_validate.m).
%
% INPUTS:
%   J   [3nV x 1] ambient current;  S surface (.Vertices, .Faces, .VertNormals, .Hemi)
%   op  (optional) a precomputed rheome.detect.operator(S), reused across frames (see rheome.detect.operator).
%
% OUTPUT (struct v):
%   .pos .nVortex .hemi .chi  (always) + .size .chirality .strength (as requested)
%
% See also: rheome.detect.criticalPoints, rheome.differential.curl, rheome.operators.connection_laplacian
%
% Author: Diellor Basha, 2026

    ALL = {'position','size','chirality','strength'};
    if nargin < 3 || isempty(props), props = ALL; end
    if ischar(props) || isstring(props), props = cellstr(props); end
    props = lower(props(:)');
    if any(strcmp(props,'all')), props = ALL; end
    props = union(props, {'position'});                  % position is the anchor -- always kept

    if nargin < 4, op = []; end
    cp = rheome.detect.criticalPoints(J, S, 'vortex', op);      % Dirac-connection detector, vortices only
    nc = numel(cp.charge);
    V  = S.Vertices;  J3 = [J(1:3:end) J(2:3:end) J(3:3:end)];
    nrm = i_normals(S);                                  % OUTWARD vertex normals
    Jt = J3 - sum(J3.*nrm,2).*nrm;  speed = vecnorm(Jt,2,2);    % tangential flow speed

    v = struct('pos', cp.pos, 'nVortex', nc, 'hemi', cp.hemi, 'chi', cp.chi);

    if any(strcmp(props,'size'))
        sz = nan(nc,1);
        for k = 1:nc
            d = vecnorm(V - cp.pos(k,:), 2, 2) * 1000;   % mm from the core (local geodesic ~ chord)
            sz(k) = i_peakradius(d, speed);              % radius (mm) where tangential speed peaks
        end
        v.size = sz;                                     % mm
    end
    if any(strcmp(props,'chirality'))
        ch = zeros(nc,1);
        for k = 1:nc
            [~,vc] = min(vecnorm(V - cp.pos(k,:), 2, 2));      % nearest vertex -> outward normal
            d = vecnorm(V - cp.pos(k,:), 2, 2) * 1000;         % mm from core
            ring = d > 3 & d < 40;                             % annulus in the swirl
            r = V(ring,:) - cp.pos(k,:);                       % displacement from core
            L = cross(r, J3(ring,:), 2);                       % r x v (rotational moment)
            ch(k) = sign(sum(L * nrm(vc,:)'));                 % (r x v).n_outward : +1 CCW, -1 CW (from outside)
        end
        v.chirality = ch;
    end
    if any(strcmp(props,'strength'))
        v.strength = cp.circulation;                     % |curl| at the core
    end
end

% ----- local helpers -----
function n = i_normals(S)
    V = S.Vertices;  nV = size(V,1);
    if isfield(S,'VertNormals') && ~isempty(S.VertNormals) && size(S.VertNormals,1)==nV
        n = S.VertNormals ./ max(vecnorm(S.VertNormals,2,2), eps);
    else
        F = double(S.Faces);  fn = cross(V(F(:,2),:)-V(F(:,1),:), V(F(:,3),:)-V(F(:,1),:));
        n = [accumarray(F(:),repmat(fn(:,1),3,1),[nV 1]), ...
             accumarray(F(:),repmat(fn(:,2),3,1),[nV 1]), ...
             accumarray(F(:),repmat(fn(:,3),3,1),[nV 1])];
        n = n ./ max(vecnorm(n,2,2), eps);
    end
end

function r = i_peakradius(d, s, win)
    if nargin < 3, win = 60; end                         % search window (mm) around the core
    edges = (0:2:win)';  m = d > 0 & d <= win;           % 2 mm rings
    b = discretize(d(m), edges);  sm = s(m);
    keep = ~isnan(b);
    mu = accumarray(b(keep), sm(keep), [numel(edges)-1 1], @mean, 0);
    ctr = edges(1:end-1) + 1;                            % ring centers (mm)
    mu = smoothdata(mu, 'movmean', 3);                   % light smoothing over rings
    [~, ip] = max(mu);  r = ctr(ip);
end

% Author: Diellor Basha, 2026
