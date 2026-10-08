function out = sensorbank(G, S, lbo, opts)
% FORWARD.SENSORBANK  A graph-wavelet bank fitted to the sensor patterns of a leadfield.
%
%   out = rheome.forward.sensorbank(G, S, lbo)
%   out = rheome.forward.sensorbank(G, S, lbo, struct('NAtoms',10,'Candidates',200))
%
% ⭐⭐ THE THREE FAMILIES ARE THE HELMHOLTZ PARTS, AND THEY SPAN THE AMBIENT 3-VECTOR. A leadfield
% row is an ambient field, three components per vertex, and at each vertex those three split into
% two tangential and one normal. Applied to a SCALAR graph wavelet w, the triple
%       grad(w)        n x grad(w)        w*n
% supplies exactly those three, so the bank is complete for smooth fields by construction. ⭐ This
% is what the Dirac operator is for elsewhere in the repo -- carrying a 3-vector on the surface --
% reached here through the Helmholtz potentials instead, which matters because a VECTOR atom has
% orientation freedom that degenerates (§41: rank 2 of 3, and a random vertex scored 0.079 against a
% fitted 0.104) while a SCALAR potential has none.
%
% ⭐⭐ THE SEARCH COSTS THREE SCALAR CORRELATIONS, NOT A VECTOR ONE. In the weak form
%       <L, grad(w)>   = -<div L, w>       <L, n x grad(w)> = <curl L, w>       <L, w*n> = <L.n, w>
% so the divergence, the curl and the normal component of the row -- three scalars on the cortex,
% each estimable directly from the standard leadfield -- score every atom by a scalar inner product.
% Only the handful of SELECTED atoms is ever materialised as a [3nV x 1] field.
%
% ⚠⚠ FIT THE FIELD, NOT THE POTENTIAL. Fitting Psi as a scalar and then differentiating loses a
% quarter: measured, a stream function fitted to 0.913 generates a solenoidal field fitted to only
% 0.664, because the gradient amplifies whatever the scalar fit left behind. This selects on the
% scalar correlation but re-solves every amplitude by least squares IN THE FIELD DOMAIN at each step,
% so the objective is always the row itself.
%
% INPUTS
%   G     [nCh x 3nV] leadfield        S  surface (.Vertices .Faces .VertNormals)
%   lbo   .Phi .Lambda .Mass           the scalar basis the wavelets are built on
%   opts  .NAtoms (10)  .Candidates (200 or a vertex vector)  .Scales (16)
%         .Sensors (all)  .Verbose (false)
%
% OUTPUT (struct out)
%   .atoms  table: sensor, iter, family, vertex, scale, coef, veField
%   .ve .curlCorr .dirCorr  [nCh x 1] against the standard leadfield row
%   .Bhat   [nCh x 3nV] the reconstructed leadfield, if opts.Reconstruct
%
% See also: rheome.differential.helmholtz, rheome.differential.curl, rheome.differential.divergence, rheome.forward.leadfield
%
% Author: Diellor Basha, 2026

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'NAtoms')     || isempty(opts.NAtoms),     opts.NAtoms = 10;     end
    if ~isfield(opts,'Candidates') || isempty(opts.Candidates), opts.Candidates = 200; end
    if ~isfield(opts,'Verbose')    || isempty(opts.Verbose),    opts.Verbose = false;  end
    if ~isfield(opts,'Reconstruct')|| isempty(opts.Reconstruct),opts.Reconstruct = false; end

    V = double(S.Vertices);  F = double(S.Faces);  nV = size(V,1);
    fg = rheome.operators.face_gradient(V, F);
    nrm = S.VertNormals ./ max(vecnorm(S.VertNormals,2,2), eps);
    Ph = lbo.Phi;  Lb = double(lbo.Lambda(:));  Ms = lbo.Mass;

    cand = opts.Candidates;
    if isscalar(cand)
        nC = round(cand);  cand = zeros(1,nC);  cand(1) = round(nV/2);
        dm = pdist2(V, V(cand(1),:));
        for i = 2:nC
            [~, cand(i)] = max(dm);  dm = min(dm, pdist2(V, V(cand(i),:)));
        end
    end
    cand = cand(:)';  nC = numel(cand);
    if ~isfield(opts,'Scales') || isempty(opts.Scales)
        lp = Lb(Lb > 1e-8*max(Lb));
        opts.Scales = 1./logspace(log10(min(lp)), log10(max(lp)), 16);
    end
    nS = numel(opts.Scales);

    % the SCALAR dictionary, built once: a graph wavelet per (candidate, scale)
    W = zeros(nV, nC*nS);  meta = zeros(nC*nS, 2);  k = 0;
    for j = 1:nC
        e = zeros(nV,1);  e(cand(j)) = 1;  ce = Ph'*(Ms*e);
        for s = 1:nS
            k = k + 1;
            gb = rheome.filters.mexhat(Lb, opts.Scales(s));
            gb = gb(:)/max(max(abs(gb)), realmin);
            W(:,k) = Ph*(gb.*ce);
            meta(k,:) = [cand(j) opts.Scales(s)];
        end
    end
    W = W - mean(W,1);
    W = W ./ max(vecnorm(W), realmin);

    gradf = @(p) fg.W*[fg.Gx*p, fg.Gy*p, fg.Gz*p];
    fam   = ["grad","rot","norm"];
    sens  = 1:size(G,1);
    if isfield(opts,'Sensors') && ~isempty(opts.Sensors), sens = opts.Sensors(:)'; end
    K = round(opts.NAtoms);

    ve = nan(size(G,1),1);  cc = ve;  dc = ve;
    rows = [];
    if opts.Reconstruct, out.Bhat = zeros(size(G)); end

    for c = sens
        L = G(c,:)';  nL = norm(L);
        if nL == 0, continue; end
        L = L/nL;
        R = L;  A = [];  pick = zeros(K,2);
        for it = 1:K
            % three scalars of the RESIDUAL drive the search
            dv = rheome.differential.divergence(R, S, fg);
            cu = rheome.differential.curl(R, S, fg);
            nu = sum(reshape(R,3,[])'.*nrm, 2);
            % ⚠⚠ THE RAW SCALAR CORRELATION IS NOT A MATCHED FILTER. <R,a> is the same quantity
            %   for all three families by the weak-form identities, but |a| is not: a gradient
            %   atom's norm grows with its scale while a normal atom's is one. Ranking on the
            %   numerator alone over-picks the rotor family -- measured 0.91 of selections against
            %   its 0.16 share of a row's energy -- and NEVER picks the normal family, which is
            %   0.34 of it. So the raw scores only SHORTLIST, and the choice is made on exact
            %   field-domain matched filters over the shortlist.
            sc = [abs(W'*dv(:)), abs(W'*cu(:)), abs(W'*nu(:))];
            nSh = 5;
            best = -inf;  a = [];  ia = 0;  fa = 0;
            for ff = 1:3
                [~, ord] = sort(sc(:,ff), 'descend');
                for jj = ord(1:min(nSh, numel(ord)))'
                    switch ff
                        case 1, av = reshape(gradf(W(:,jj))', [], 1);
                        case 2, av = reshape(cross(nrm, gradf(W(:,jj)), 2)', [], 1);
                        case 3, av = reshape((W(:,jj).*nrm)', [], 1);
                    end
                    na = norm(av);
                    if na == 0, continue; end
                    r = abs(av'*R)/na;                 % ⭐ the exact matched filter
                    if r > best, best = r; a = av; ia = jj; fa = ff; end
                end
            end
            if isempty(a) || norm(a) == 0, break; end
            A = [A a/norm(a)]; %#ok<AGROW>
            pick(it,:) = [ia fa];
            Q = orth(A);
            R = L - Q*(Q'*L);            % ⭐ re-solve in the FIELD domain every step
        end
        Q = orth(A);  Lhat = Q*(Q'*L);
        ve(c) = 1 - norm(R)^2/norm(L)^2;
        cL = rheome.differential.curl(L, S, fg);  cA = rheome.differential.curl(Lhat, S, fg);
        cc(c) = corr(cL(:), cA(:));
        L3 = reshape(L,3,[])';  A3 = reshape(Lhat,3,[])';  wgt = vecnorm(L3,2,2);
        u = L3./max(vecnorm(L3,2,2),realmin);  z = A3./max(vecnorm(A3,2,2),realmin);
        dc(c) = sum(wgt.*sum(u.*z,2))/sum(wgt);
        if opts.Reconstruct, out.Bhat(c,:) = (Lhat*nL)'; end
        cf = Q'*L;
        for it = 1:size(pick,1)
            if pick(it,1) == 0, continue; end
            rows = [rows; c it double(pick(it,2)) meta(pick(it,1),1) meta(pick(it,1),2) ...
                    cf(min(it,numel(cf)))]; %#ok<AGROW>
        end
        if opts.Verbose, fprintf('  sensor %d: ve %.3f\n', c, ve(c)); end
    end

    out.ve = ve;  out.curlCorr = cc;  out.dirCorr = dc;
    out.atoms = table(rows(:,1), rows(:,2), fam(rows(:,3))', rows(:,4), rows(:,5), rows(:,6), ...
        'VariableNames', {'sensor','iter','family','vertex','scale','coef'});
    out.candidates = cand;  out.scales = opts.Scales;
end

% Author: Diellor Basha, 2026
