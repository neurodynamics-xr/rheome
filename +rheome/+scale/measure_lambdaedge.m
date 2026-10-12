function [T, X] = measure_lambdaedge(name, S, opts)
% SCALE.MEASURE_LAMBDAEDGE  The shortest wavelength at which upstream readout v2 reads travel along G1 the right way, on one head.
%
%   [T, X] = rheome.scale.measure_lambdaedge(name)
%   [T, X] = rheome.scale.measure_lambdaedge(name, S, Gradient='/path/g1_fsaverage.mat')
%   [T, X] = rheome.scale.measure_lambdaedge(name, S, Lambda=[150 200 250 300 400 500 700 1000])
%
% Why: readout v2 (rheome.scale.measure_upstream, Readout="v2") reads planted travel along the principal
% gradient right down to ~200 mm in some heads and only from ~500 mm in others; the short-wavelength edge is a
% property of the head (its folding and sensors), not of a recording. A result read through v2 states the
% wavelengths it covers, so each head's edge is measured, and from the head alone.
%
% What: plane waves exp(-i s 2 pi x / lambda) of normal current along the participant's G1 axis x (the
% area-weighted linear fit of the signed G1 map on the vertex positions), s = +1 (up the gradient) and -1 (down),
% noiseless, through the participant's gain. Each is read by v2 exactly as measure_upstream reads tiles (its
% geometry: depth-3 tiles, fwdLS, G1 per tile, own-tile mask) and by the tile means of the planted current
% itself (truth, no sensors), with rheome.flow.upstreamshare over four phase frames, the hemispheres pooled by
% their weight.
%   right   a lambda reads right when share_v2 > 0.5 for s = +1 and < 0.5 for s = -1.
%   EDGE    edge_mm = the shortest grid lambda from which v2 reads right at EVERY grid lambda up to the
%           largest (1000 mm by default); NaN when v2 does not read right at the largest.
%
% NO RECORDING IS READ: the inputs are the cortex, channels, gain, noise covariance (for the kernel and its
% own-tile mask) and the G1 map. Nothing here is an outcome: there is no U, U_surrogate or U* of any recording.
%
% Returns T (rheome.scale.rows, analysis "lambdaedge"): edge_mm, readout_version (2), samples_read (samples in
% the subject's cached recording: 0 under rheome.scale.run Import="head"), kept_weight_L/R, g1_tilecheck_L/R
% and one lambda_grid_mm row per grid lambda; and X, one row per (lambda_mm, dir): share_truth, share_v2,
% right_truth, right_v2, v2_phase_err_med_deg (median over tiles of |angle(q conj(truth))|) and in_grid (false
% for the Reported lambdas, which do not enter the edge).
%
% ⚠ The tiles are ~125 mm across: v2 aliases at lambda ~100 mm in most heads, and below ~240 mm any tile
%   readout can read backwards. The plant has v2's own gain and normals (the inverse crime), so the edge is a
%   best case for that head; it says nothing about noise.
%
% See also: rheome.scale.measure_upstream, rheome.flow.upstreamshare, rheome.scale.run
%
% Author: Diellor Basha, 2026

    arguments
        name (1,:) char
        S = []
        opts.Gradient char = ''              % default: RHEOME_G1
        opts.SignCheck (1,1) logical = true
        opts.OwnMin (1,1) double = 0.25     % v1/v2: 0.25; lower it only on a toy head (tests)
        opts.Lambda (1,:) double {mustBePositive} = [150 200 250 300 400 500 700 1000]   % mm: the edge's grid
        opts.Reported (1,:) double {mustBePositive} = [1500 2000]                       % mm: reported, not in the edge
    end
    if isempty(S), S = rheome.scale.sensors(name); end
    H = rheome.scale.measure_upstream(name, S, Readout="v2", Geometry=true, Gradient=opts.Gradient, ...
        SignCheck=opts.SignCheck, OwnMin=opts.OwnMin);
    lam = sort(opts.Lambda);
    V = double(S.Sf.Vertices);  N = double(S.Sf.VertNormals);  N = N ./ vecnorm(N, 2, 2);
    Gn = N(:,1)' .* S.G(:, 1:3:end) + N(:,2)' .* S.G(:, 2:3:end) + N(:,3)' .* S.G(:, 3:3:end);   % normal-current gain
    w = zeros(S.nV, 1);  for h = ["L" "R"], w(H.(h).gv) = H.(h).a; end
    gv = [H.L.gv; H.R.gv];  gg = H.sgn * [H.L.gv_g; H.R.gv_g];  ok = isfinite(gg);  ww = sqrt(w(gv(ok)));
    b = ([ones(nnz(ok), 1) V(gv(ok), :)] .* ww) \ (gg(ok) .* ww);  x = V * (b(2:4) / norm(b(2:4)));   % G1 axis, m

    fr = reshape(exp(1i * 2*pi * (0:3) / 4), 1, 1, []);
    X = table();
    for l = [lam sort(opts.Reported)]
        for s = [1 -1]
            c = exp(-1i * s * 2*pi * x / (l / 1000));  q = H.Pinv * (Gn * (w .* c));
            num = [0 0];  den = [0 0];  pe = [];  off = 0;
            for h = ["L" "R"]                                            % H.Pinv's rows: L tiles, then R tiles
                Hh = H.(h);  n = numel(Hh.g);  tru = Hh.Pa' * c(Hh.gv);  qh = q(off + (1:n));  off = off + n;
                o = [rheome.flow.upstreamshare(tru .* fr, Hh.g, Hh.A, Mask=Hh.mask), ...
                     rheome.flow.upstreamshare(qh .* fr, Hh.g, Hh.A, Mask=Hh.mask)];
                num = num + [o.U] .* [o.weight];  den = den + [o.weight];
                pe = [pe; abs(angle(qh .* conj(tru)))]; %#ok<AGROW>
            end
            u = num ./ den;
            X = [X; table(l, s, u(1), u(2), (u(1) - 0.5) * s > 0, (u(2) - 0.5) * s > 0, median(rad2deg(pe)), ...
                ismember(l, lam), 'VariableNames', {'lambda_mm', 'dir', 'share_truth', 'share_v2', 'right_truth', ...
                'right_v2', 'v2_phase_err_med_deg', 'in_grid'})]; %#ok<AGROW>
        end
    end
    edge = i_edge(lam, arrayfun(@(l) all(X.right_v2(X.lambda_mm == l)), lam));

    st = rheome.load.study(name);  nRead = size(st.rec.F, 2);  clear st
    T = [rheome.scale.rows("lambdaedge", ["edge_mm" "readout_version" "samples_read" "kept_weight_L" "kept_weight_R" ...
             "g1_tilecheck_L" "g1_tilecheck_R"], [edge 2 nRead H.L.keptWeight H.R.keptWeight H.L.tilecheck H.R.tilecheck], ...
             ["mm" "version" "samples" "fraction" "fraction" "flag" "flag"]);
         rheome.scale.rows("lambdaedge", repmat("lambda_grid_mm", numel(lam), 1), lam, "mm")];
    fprintf('[lambdaedge %s] v2 edge %g mm (right at %s mm)\n', name, edge, strjoin(string(lam(arrayfun(@(l) ...
        all(X.right_v2(X.lambda_mm == l)), lam))), ' '));
end

function e = i_edge(lam, right)
% the shortest grid lambda from which every grid lambda up to the largest reads right; NaN if the largest does not
    k = find(~right, 1, 'last');
    if isempty(k), e = lam(1); elseif k == numel(lam), e = NaN; else, e = lam(k + 1); end
end

% Author: Diellor Basha, 2026
