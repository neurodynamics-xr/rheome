function out = sensoratom(G, dbasis, opts)
% FORWARD.SENSORATOM  Fit each sensor's sensitivity as a GRAPH WAVELET on the cortex.
%
%   out = rheome.forward.sensoratom(G, dbasis)
%   out = rheome.forward.sensoratom(G, dbasis, struct('Surface',S,'Candidates',300,'NAtoms',3))
%
% ⭐⭐ A DICTIONARY DERIVED FROM THE INSTRUMENT INSTEAD OF ASSUMED. rheome.flow.vortexbank is a bank of
% atoms chosen by hand, and §38 measured the consequence: the parametric inverse beats MNE only when
% the source happens to match that choice, and is 44 mm wrong when it does not. The atoms fitted here
% are whatever the leadfield actually is, so the dictionary carries no prior about source shape.
%
% THE FORM, which falls straight out of the mode representation. A graph wavelet centred at vertex v
% with spectral profile g has mode coefficients g(lambda_k)*d_v(k), where d_v is the delta at v
% expanded in the basis. So a sensor row's coefficients y are fitted by choosing v, choosing g from a
% kernel family, and solving for the three orientation weights in closed form -- the atom at a vertex
% spans a 3-D subspace, one per Cartesian component of the delta.
%
% ⭐ MEASURED, ON REAL MEG (270 sensors, 400 Dirac modes, docs section 40):
%     one atom     var explained  0.688     three atoms  0.789 - 0.827
%     a VORTEX atom, for contrast   0.087   five vortices 0.104
%   So a sensor IS a wavelet over the eigenmodes, and it is not a vortex: the same row that a vortex
%   explains 9% of, a plain graph wavelet explains 69% of.
%
% ⚠⚠ THE CENTRE IS NOT UNDER THE SENSOR. Median distance from the fitted centre to the cortical
%   vertex directly beneath the sensor is 56 mm at a 12 mm candidate spacing. That is physics rather
%   than misfit: MEG is blind to radial sources and the cortex directly beneath a sensor is radial
%   there, so peak sensitivity sits on the flanks. Do not read `vertex` as "the sensor's location".
%
% ⚠ Refining the candidate grid from 24 mm to 12 mm moves one atom from 0.647 to 0.688 but moves
%   THREE atoms DOWN, 0.827 to 0.789 -- a finer grid lets the greedy triple pick three near-identical
%   centres. Fit the number of atoms you intend to use.
%
% INPUTS
%   G       [nCh x 3nV] vertex leadfield (rheome.forward.leadfield)
%   dbasis  .Phi [4nV x nModes], .Lambda, .nVert, .nModes -- as rheome.load.dirac restricted to a hemisphere
%   opts    .Surface    the hemisphere surface, needed to pick candidates by farthest point
%           .Candidates count (300) or an explicit vertex vector
%           .Scales     t values for the kernel (25 logarithmic, spanning the spectrum)
%           .Family     "mexhat" (default) | "heat"
%           .NAtoms     atoms per sensor (1)
%           .Verbose    false
%
% OUTPUT (struct out)
%   .table   one row per sensor: sensor, vertex, tau, varExplained
%   .vertex .tau .varExplained  [nCh x NAtoms] / [nCh x 1]
%   .candidates  the vertices searched      .family
%
% See also: rheome.forward.diracgain, rheome.forward.leadfield, rheome.flow.vortexbank, rheome.filters.mexhat, rheome.filters.heat
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'Candidates') || isempty(opts.Candidates), opts.Candidates = 300;  end
    if ~isfield(opts,'Family')     || isempty(opts.Family),     opts.Family     = "mexhat"; end
    if ~isfield(opts,'NAtoms')     || isempty(opts.NAtoms),     opts.NAtoms     = 1;    end
    if ~isfield(opts,'Verbose')    || isempty(opts.Verbose),    opts.Verbose    = false; end
    assert(any(strcmpi(opts.Family, ["mexhat","heat"])), ...
        'rheome.forward.sensoratom: Family must be "mexhat" or "heat".');

    lam  = double(dbasis.Lambda(:));
    lmax = max(lam);
    assert(lmax > 0, 'rheome.forward.sensoratom: the basis has no positive eigenvalue.');
    [Gm, Psi3] = rheome.forward.diracgain(G, dbasis);
    nCh = size(Gm,1);  nV = dbasis.nVert;

    if ~isfield(opts,'Scales') || isempty(opts.Scales)
        opts.Scales = logspace(0, 8, 25)/lmax;
    end

    cand = opts.Candidates;
    if isscalar(cand)
        assert(isfield(opts,'Surface') && ~isempty(opts.Surface), ...
            'rheome.forward.sensoratom: a Candidates COUNT needs opts.Surface to sample from.');
        V = double(opts.Surface.Vertices);
        nC = round(cand);  cand = zeros(1,nC);  cand(1) = round(nV/2);
        dm = pdist2(V, V(cand(1),:));
        for i = 2:nC
            [~, cand(i)] = max(dm);
            dm = min(dm, pdist2(V, V(cand(i),:)));
        end
    end
    cand = cand(:)';

    % every (candidate, scale) subspace, precomputed once and shared by all sensors
    nCand = numel(cand);  nS = numel(opts.Scales);
    Qs = cell(nCand, nS);
    for j = 1:nCand
        D3 = Psi3(3*(cand(j)-1)+(1:3), :)';
        for s = 1:nS
            if strcmpi(opts.Family,"mexhat"), gk = rheome.filters.mexhat(lam, opts.Scales(s));
            else,                             gk = rheome.filters.heat(lam, opts.Scales(s), lmax); end
            % ⚠⚠ NORMALISE THE KERNEL OR THE ATOM'S RANK COLLAPSES. At the scales this fit
            %   selects, the mexhat peaks below the smallest eigenvalue and its values run to
            %   1e-11 and below. orth() uses a tolerance relative to the largest singular value
            %   of the whole matrix, so the three columns -- which differ only in the delta's
            %   orientation -- were being discarded: measured rank 1 of 3 on real data, and the
            %   score was then computed in a one-dimensional subspace. Scaling is a no-op
            %   mathematically, since the coefficients absorb it, and decisive numerically.
            gk = gk(:) / max(max(abs(gk)), realmin);
            Q = orth(gk.*D3);
            if ~isempty(Q), Qs{j,s} = Q; end
        end
        if opts.Verbose && mod(j,50)==0
            fprintf('  sensoratom: %d/%d candidates\n', j, nCand);
        end
    end

    K = round(opts.NAtoms);
    vtx = zeros(nCh, K);  tau = zeros(nCh, K);  ve = zeros(nCh, 1);
    for c = 1:nCh
        y = Gm(c,:)';  ny = norm(y);
        if ny == 0, continue; end
        y = y/ny;
        sc = -inf(nCand, 1);  st = ones(nCand, 1);
        for j = 1:nCand
            for s = 1:nS
                if isempty(Qs{j,s}), continue; end
                r = norm(Qs{j,s}'*y)^2;
                if r > sc(j), sc(j) = r; st(j) = s; end
            end
        end
        [~, ord] = sort(sc, 'descend');
        pick = ord(1:min(K, nCand));
        Qall = [];
        for i = 1:numel(pick)
            Qall = [Qall Qs{pick(i), st(pick(i))}]; %#ok<AGROW>
            vtx(c,i) = cand(pick(i));  tau(c,i) = opts.Scales(st(pick(i)));
        end
        Qo = orth(Qall);
        ve(c) = norm(Qo'*y)^2;
    end

    out.vertex = vtx;  out.tau = tau;  out.varExplained = ve;
    out.candidates = cand;  out.family = string(opts.Family);  out.scales = opts.Scales;
    out.table = table((1:nCh)', vtx(:,1), tau(:,1), ve, ...
        'VariableNames', {'sensor','vertex','tau','varExplained'});
end

% Author: Diellor Basha, 2026
