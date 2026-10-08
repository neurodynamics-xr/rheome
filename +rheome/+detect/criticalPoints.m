function cp = criticalPoints(J, S, types, op)
% DETECT.CRITICALPOINTS  Classified critical points of an ambient source vector field.
%
%   cp = rheome.detect.criticalPoints(J, S)               % all types
%   cp = rheome.detect.criticalPoints(J, S, types)        % only the requested type(s)
%   cp = rheome.detect.criticalPoints(J, S, types, op)    % reuse a precomputed rheome.detect.operator(S)
%
% The Dirac-connection route, chained end to end and per hemisphere:
%   1. rheome.operators.connection_laplacian  -> Levi-Civita connection + vertex tangent frames
%   2. z = (J.e1) + i(J.e2)             -> the complex tangent field
%   3. winding number per triangle      -> integer TOPOLOGICAL CHARGE (Poincare-Hopf exact)
%   4. CLASSIFY each charge with the local divergence & curl (rheome.differential.*): the winding is
%      +1 for a VORTEX and a SOURCE/SINK alike, so only div/curl tells them apart (validated
%      on analytic fields in detect_validate.m).
%
% Type (from charge + local div/curl):
%   charge +1, |curl| >  |div|         -> 'vortex'   (rotation: centre / spiral)
%   charge +1, |curl| <= |div|, div>0  -> 'source'
%   charge +1, |curl| <= |div|, div<0  -> 'sink'
%   charge -1                          -> 'saddle'
%
% INPUTS:
%   J      [3nV x 1] ambient current; S surface (.Vertices, .Faces, .Hemi to split hemispheres)
%   types  which to return: a char ('vortex') or cellstr ({'vortex','saddle'}); also 'all'
%          (default) and 'node' (= source + sink). Case-insensitive.
% ⚠⚠ THIS LOCALISES SINGULARITIES TO FACES, so it cannot confirm a singularity that was PRESCRIBED AT A
% VERTEX. Direction-field design (rheome.flow.directionfield, Crane-Desbrun-Schroeder) plants index k_v at
% vertices; such a singularity has zero winding on every individual triangle around that vertex while
% carrying index k on the vertex LINK. Measured: with +1 prescribed at two vertices, the vertex-link
% index reads +1 and 0 while this function puts 42 nonzero faces nowhere near either. ⭐ Validate a
% DESIGNED field by its holonomy, which is exact; use this for a MEASURED field, which is what it is for.
% ⚠ It also cannot report |charge| > 1 -- the selection is `find(abs(q) == 1)` -- so a higher-order
% defect comes back as +-1.
%
%   op     (optional) a precomputed rheome.detect.operator(S): the whole-surface face gradient + the
%          per-hemisphere connection Laplacian, built once and reused. Omit -> built internally
%          (identical result). Pass it in per-frame loops to skip rebuilding the operators.
%
% OUTPUT (struct cp, arrays [nc x 1] sorted by DESCENDING strength = max(|curl|,|div|)):
%   .pos [nc x 3]  .charge (+1/-1)  .type {cellstr}  .circulation (|curl|)  .divergence (signed)
%   .strength (max(|curl|,|div|))  .hemi  .chi [1 x nHemi]  (Poincare-Hopf check = Euler chi)
%
% See also: rheome.operators.connection_laplacian, rheome.differential.divergence, rheome.differential.curl, detect_validate
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(types), types = {'vortex','source','sink','saddle'}; end
    if ischar(types) || isstring(types), types = cellstr(types); end
    types = lower(types(:)');
    if any(strcmp(types,'all')),  types = {'vortex','source','sink','saddle'}; end
    if any(strcmp(types,'node')), types = [types, {'source','sink'}]; end

    if nargin < 4 || isempty(op), op = rheome.detect.operator(S); end          % build once if not supplied
    if op.nV ~= size(J,1)/3
        error('detect:criticalPoints:opMismatch', 'op is for %d vertices but J has %d.', op.nV, size(J,1)/3);
    end
    nH    = op.nH;
    J3    = [J(1:3:end), J(2:3:end), J(3:3:end)];
    divv  = rheome.differential.divergence(J, S, op.fg);   % [nV] whole-surface divergence (source +, sink -)
    curlv = rheome.differential.curl(J, S, op.fg);         % [nV] whole-surface vorticity

    pos=zeros(0,3); charge=zeros(0,1); circ=zeros(0,1); dvg=zeros(0,1); hemi=zeros(0,1); ty={};
    chi = zeros(1, nH);
    for hh = 1:nH
        Sh = op.Sh{hh};  gv = op.gv{hh};  C = op.C{hh};    % precomputed hemisphere + connection
        Jl = J3(gv, :);
        z  = sum(Jl.*C.e1, 2) + 1i*sum(Jl.*C.e2, 2);

        F = double(Sh.Faces);  a=F(:,1); b=F(:,2); c=F(:,3);  wr = @(x) mod(x+pi,2*pi)-pi;
        idx = ( wr(angle(z(b))-angle(z(a))-angle(i_edge(C.Rt,a,b))) ...
              + wr(angle(z(c))-angle(z(b))-angle(i_edge(C.Rt,b,c))) ...
              + wr(angle(z(a))-angle(z(c))-angle(i_edge(C.Rt,c,a))) ) / (2*pi);
        q = round(idx);  chi(hh) = sum(q);

        sel = find(abs(q) == 1);
        gf  = [gv(a(sel)), gv(b(sel)), gv(c(sel))];
        dC  = mean(divv(gf), 2);   cC = mean(curlv(gf), 2);
        cen = (Sh.Vertices(a(sel),:) + Sh.Vertices(b(sel),:) + Sh.Vertices(c(sel),:)) / 3;
        for k = 1:numel(sel)
            if     q(sel(k)) == -1,          tk = 'saddle';
            elseif abs(cC(k)) >  abs(dC(k)), tk = 'vortex';
            elseif dC(k) > 0,                tk = 'source';
            else,                            tk = 'sink';
            end
            ty{end+1,1} = tk; %#ok<AGROW>
        end
        pos=[pos;cen]; charge=[charge;q(sel)]; circ=[circ;abs(cC)]; dvg=[dvg;dC]; hemi=[hemi;hh*ones(numel(sel),1)]; %#ok<AGROW>
    end

    keep = ismember(ty, types);                             % filter to the requested type(s)
    pos=pos(keep,:); charge=charge(keep); circ=circ(keep); dvg=dvg(keep); hemi=hemi(keep); ty=ty(keep);
    strength = max(circ, abs(dvg));                         % type-agnostic saliency
    [strength, o] = sort(strength, 'descend');
    cp = struct('pos',pos(o,:), 'charge',charge(o), 'type',{ty(o)}, 'circulation',circ(o), ...
                'divergence',dvg(o), 'strength',strength, 'hemi',hemi(o), 'chi',chi);
end

% ----- helper -----
function v = i_edge(Rt, r, c)
    v = full(Rt(sub2ind(size(Rt), r, c)));
end

% Author: Diellor Basha, 2026
