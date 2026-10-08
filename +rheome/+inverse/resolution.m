function out = resolution(K, Gain, basis, opts)
% INVERSE.RESOLUTION  What an inverse can actually resolve on the cortex, in millimetres.
%
%   out = rheome.inverse.resolution(Res.ImagingKernel, Gain, B.L.lbo, Vertices=B.L.S.Vertices, ...
%                            Faces=B.L.S.Faces, GlobalIdx=B.L.gv)
%   out = rheome.inverse.resolution(K, Gain, lbo, ..., NumSeeds=300, Modes=1:1000)
%
% Reduces the resolution matrix R = K*Gain to LENGTHS. Two families of number, because they
% answer two different questions and they do NOT agree:
%
%   LOCALISATION (.r50 .sd .ple .offHemi)   how far a point source's estimate spreads.
%       Bounds "are these two blobs separate" and "how big is this vortex".
%   SPATIAL FREQUENCY (.mtf .mtfBand)       the gain applied to a periodic cortical pattern.
%       Bounds "can a pattern of this wavelength contribute at all".
%
% ⚠ MEASURED ON A REFERENCE SUBJECT, 270 CHANNELS, SnrFixed = 3: r50 = 52 mm (51.6, IQR 42.5-76.1;
% equivalent wavelength 262 mm) against an MTF half-maximum of 99 mm -- a factor of 2.6.
% ⚠ .mtfHalf itself is unstable -- bimodal across resting cohorts (most subjects 26-44 mm, a
% minority 57-111 mm, the reference subject among the minority) -- so prefer .mtfCentroid for any
% comparison. A system can localise badly
% and still respond to a fine GLOBAL grating, because the grating is a coherent sum over the
% whole surface while a point source is not. Quote the one that matches the claim being made.
% Anything about the size, count or separateness of spatial features is a LOCALISATION claim.
%
% ⭐ trace(R) IS THE EFFECTIVE RANK, and it is the same number the regularisation sets.
% R = VL*diag(g.*s)*VL' with g.*s = Lambda*s^2/(Lambda*s^2+1), so the Wiener gains of
% rheome.inverse.mne stage 3 ARE the eigenvalues of the resolution matrix: .dof = sum of them.
% On this array .dof is 48.7 at SnrFixed 2.2, 57.5 at 3 and 113.4 at 15 -- of 270 channels.
%
% ⭐ THE r50 -> WAVELENGTH CONVERSION IS MEASURED, NOT ASSUMED (Calibrate=true, the default).
% Heat-kernel blobs of known aperture are planted on the same mesh with the same seeds and the
% same statistic, giving r50 = 0.196*wavelength at R^2 = 0.9993 over 40-350 mm on this cortex.
% A conversion argued from a Gaussian instead would be a definition; this is a measurement.
%
% ⚠ THE CONVERSION SATURATES WHEN THE APERTURE APPROACHES THE PATCH ITSELF. Geodesic distance
% on a closed surface stops growing at its own diameter, so r50 stops growing too and the line
% bends over. Measured: on the reference hemisphere (max geodesic 314 mm) apertures of 40-350 mm
% give R^2 = 0.999, while on a 70 mm sphere the same range gives 0.80. Keep Apertures well
% under the geodesic diameter and check .r50R2 before trusting .lambdaLoc.
%
% ⚠ THE MTF LIVES ON A SUBSPACE, AND IT IS A SMALL ONE. The normally oriented modes n*phi_k of
% ONE hemisphere are orthonormal in the mass metric (|n| = 1) but span nV of 3nV dimensions on
% half the cortex, and they are WORSE observed than that share suggests: their trace is 5.8 of
% R's 57.5 at SnrFixed 3, about 10% where the naive share is 17%. A smooth normally oriented
% pattern flips direction across every gyrus and largely cancels at the sensors. So .mtf is the
% response to normally oriented cortical patterns, which is the physiological case, not to
% every current R can see.
%
% ⚠⚠ IT IS NOT MONOTONE, AND .mtfHalf CAN THEREFORE JUMP. On this cortex the 31-62 mm octave
% is recovered slightly BETTER than 62-125 mm (0.49 against 0.45 of the peak), because the
% folding itself has a wavelength there: a normally oriented pattern at gyral scale is the one
% that does not cancel. So .mtfHalf -- the first crossing below half going finer -- sits just
% above that dip, and any change that lifts the 62-125 octave past 0.5 sends the crossing
% clean past it to the next octave. MEASURED: raising SnrFixed from 2.2 to 51.3 moves that
% octave from 0.46 to 0.54, a change of 0.08, and .mtfHalf reports 96 mm then 33 mm -- a
% threefold "improvement" that is entirely an artefact of reading a threshold off a bumpy
% curve. ⭐ QUOTE .mtfCentroid INSTEAD when comparing curves; it is the gain-weighted centroid
% in log wavelength and cannot jump. Read .mtfBand as a curve either way.
%
% ⚠ THE SHAPE OF .mtf BARELY MOVES WITH SnrFixed while .dof more than doubles. Regularisation
% sets HOW MUCH is recovered, not WHICH wavelengths: the normalised band curve changes by under
% 0.1 between SnrFixed 2.2 and 15. Do not expect a better SNR to buy finer spatial detail.
%
% INPUTS:
%   K       [3V x C] imaging kernel, whitener folded in (rheome.inverse.mne -> .ImagingKernel).
%           ⚠ USE 'amplitude' MEASURE. A dSPM or sLORETA kernel is already divided by a
%           per-vertex noise normaliser, so R = K*Gain is no longer the resolution matrix.
%   Gain    [C x 3V] the SAME channel rows the kernel was built from
%   basis   cached per-hemisphere lbo struct (.Phi [nH x Ks], .Lambda, .Mass) -- rheome.load.bases,
%           NOT an eigensolve
%   opts (name-value):
%     Vertices, Faces   the hemisphere's own surface (rheome.load.bases -> B.L.S)
%     Normals           [nH x 3] vertex normals; pass B.L.S.VertNormals rather than letting
%                       them be rebuilt from the faces
%     GlobalIdx         [1 x nH] this hemisphere's rows in the whole-cortex vertex numbering
%                       (B.L.gv). ⚠ REQUIRED when the kernel covers both hemispheres.
%     NumSeeds          point sources for the localisation family (300); 0 skips it
%     Modes             mode indices for the MTF ([] = all of them); 0 skips it
%     Calibrate         measure the r50 -> wavelength conversion (true)
%     Apertures         planted wavelengths for that calibration, metres
%
% OUTPUT (struct):
%   .r50 .sd .ple .offHemi [nSeeds x 1]   geodesic radius holding half the PSF power; spatial
%           dispersion sqrt(E[d^2]); peak localisation error; fraction of power on the WRONG
%           hemisphere -- all in metres except the last
%   .seeds  local vertex indices    .depth [nSeeds x 1] distance to the nearest given coil, if
%           SensorLoc was passed
%   .r50PerWavelength  the measured conversion c in r50 = c*wavelength, and .r50R2
%   .lambdaLoc         median r50 / c: the localisation floor AS A WAVELENGTH
%   .mtf [nModes x 1]  <n phi_k, R n phi_k>, the per-mode gain     .modeWavelength metres
%   .mtfBand .mtfBandEdges .mtfHalf .mtfPeak .mtfCentroid   octave-averaged curve, its edges,
%           the half-max crossing going finer (⚠ jumpy: see above), the centre of the peak
%           octave, and the gain-weighted centroid in log wavelength (⭐ the robust one)
%   .dofSubspace  sum(.mtf): how much of the rank the tested subspace carries
%
% See also: rheome.inverse.mne, rheome.flow.crosstalk, rheome.flow.sensitivity, rheome.sensors.calibrate, rheome.geom.tree
%
% Author: Diellor Basha, 2026

    arguments
        K       double
        Gain    double
        basis   (1,1) struct
        opts.Vertices   double = []
        opts.Faces      double = []
        opts.GlobalIdx  double = []
        opts.Normals    double = []
        opts.SensorLoc  double = []
        opts.NumSeeds   (1,1) double {mustBeNonnegative} = 300
        opts.Modes      double = []
        opts.Calibrate  (1,1) logical = true
        opts.Apertures  double = [40 60 90 120 180 250 350]*1e-3
        opts.Seed       (1,1) double = 7
    end
    Phi = basis.Phi;  Lam = basis.Lambda(:);  Mss = basis.Mass;
    nH  = size(Phi, 1);
    V   = size(Gain, 2) / 3;
    gv  = opts.GlobalIdx;  if isempty(gv), gv = 1:nH; end
    gv  = double(gv(:))';
    if numel(gv) ~= nH
        error('inverse:resolution:idx', 'GlobalIdx has %d entries but the basis has %d vertices.', ...
              numel(gv), nH);
    end
    if size(K,2) ~= size(Gain,1)
        error('inverse:resolution:chan', 'Kernel has %d channels, leadfield has %d rows.', ...
              size(K,2), size(Gain,1));
    end
    rows = reshape((gv-1)*3 + (1:3)', [], 1);       % this hemisphere's rows in the 3V source space
    mv   = full(sum(Mss, 2));
    out  = struct();

    % ---- geodesic metric on the hemisphere: Dijkstra with Euclidean edge weights ----
    doLoc = opts.NumSeeds > 0;
    if doLoc
        if isempty(opts.Vertices) || isempty(opts.Faces)
            error('inverse:resolution:surf', 'Vertices and Faces are needed for the localisation family.');
        end
        P = opts.Vertices;  F = opts.Faces;
        E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
        w = vecnorm(P(E(:,1),:) - P(E(:,2),:), 2, 2);
        Gg = simplify(graph(E(:,1), E(:,2), w, nH));
        rs = rng(opts.Seed);  cleanupRng = onCleanup(@() rng(rs));
        ns = min(opts.NumSeeds, nH);
        seeds = sort(randperm(nH, ns))';
        D = distances(Gg, seeds);                   % [ns x nH]

        r50 = nan(ns,1); sd = nan(ns,1); ple = nan(ns,1); offh = nan(ns,1);
        for j = 1:ns
            v  = gv(seeds(j));
            Eb = K * Gain(:, (v-1)*3 + (1:3));      % [3V x 3] column block of R
            p  = sum(reshape(sum(Eb.^2, 2), 3, []), 1)';   % orientation-free PSF power
            pL = p(gv);
            offh(j) = 1 - sum(pL)/sum(p);
            [r50(j), sd(j), ple(j)] = i_spread(pL, D(j,:)');
        end
        out.seeds = seeds;  out.r50 = r50;  out.sd = sd;  out.ple = ple;  out.offHemi = offh;
        if ~isempty(opts.SensorLoc)
            out.depth = min(pdist2(P(seeds,:), opts.SensorLoc'), [], 2);
        end

        % ---- the conversion, measured on planted blobs with the same seeds and statistic ----
        if opts.Calibrate
            wp = opts.Apertures(:);  rp = nan(numel(wp),1);
            for i = 1:numel(wp)
                gk = exp(-Lam * (wp(i)/(2*pi))^2);
                rr = nan(ns,1);
                for j = 1:ns
                    d0 = zeros(nH,1);  d0(seeds(j)) = 1;
                    x  = Phi * (gk .* (Phi' * (Mss * d0)));
                    rr(j) = i_spread(x.^2, D(j,:)');
                end
                rp(i) = median(rr);
            end
            c = (rp' * wp) / (wp' * wp);            % through the origin: r50 = c * wavelength
            out.r50PerWavelength = c;
            out.r50R2  = 1 - sum((rp - c*wp).^2) / sum((rp - mean(rp)).^2);
            out.lambdaLoc = median(r50) / c;
            out.calibration = table(wp, rp, rp./wp, 'VariableNames', {'wavelength','r50','ratio'});
        end
    end

    % ---- the MTF: the diagonal of R in the normally oriented mode family ----
    if ~isequal(opts.Modes, 0)
        md = opts.Modes;  if isempty(md), md = 1:numel(Lam); end
        md = double(md(:))';
        Nrm = i_normals(opts, nH);
        d = nan(numel(md), 1);  blk = 200;
        for a = 1:blk:numel(md)
            j = md(a:min(a+blk-1, numel(md)));
            X = zeros(3*V, numel(j));
            for c = 1:numel(j), X(rows, c) = reshape((Nrm .* Phi(:, j(c)))', [], 1); end
            Y = K * (Gain * X);
            for c = 1:numel(j)
                jn = sum(reshape(Y(rows, c), 3, [])' .* Nrm, 2);
                d(a+c-1) = sum(mv .* Phi(:, j(c)) .* jn);
            end
            clear X Y
        end
        out.mtf = d;  out.modes = md(:);
        out.modeWavelength = 2*pi ./ sqrt(max(Lam(md), eps));
        out.dofSubspace = sum(d);
        % ⚠ THE CONSTANT MODE HAS NO WAVELENGTH. Lambda(1) is ~1e-12 here, so 2*pi/sqrt(Lambda)
        % is 6e6 m and it would set the top octave edge by itself. Only the numerically zero
        % eigenvalues -- one per connected component -- are excluded from the banding, by a
        % tolerance on the eigenvalue and NOT by a bound on the wavelength: the first genuine
        % mode of this hemisphere is at 545 mm against a 378 mm equivalent diameter, which is
        % correct for a fundamental and must not be filtered out. They stay in .mtf and in
        % .dofSubspace, where they belong.
        keep = Lam(md) > max(Lam) * 1e-9;
        out.mtfBandExcluded = sum(~keep);
        [out.mtfBand, out.mtfBandEdges, out.mtfHalf, out.mtfPeak, out.mtfCentroid] = ...
            i_bands(d(keep), out.modeWavelength(keep));
    end
end

function [r50, sd, ple] = i_spread(p, d)
% The three localisation numbers from one PSF and its geodesic distances.
    p = max(p, 0);  s = sum(p);
    if s <= 0, r50 = NaN; sd = NaN; ple = NaN; return; end
    q = p / s;
    [~, im] = max(p);  ple = d(im);
    sd = sqrt(sum(q .* d.^2));
    [ds, o] = sort(d);  cq = cumsum(q(o));
    r50 = ds(find(cq >= 0.5, 1));
end

function [bnd, edges, half, peak, centroid] = i_bands(d, lam)
% Octave-average the per-mode gain, then read the half-maximum crossing going finer.
    top = 2^ceil(log2(max(lam)));
    edges = top ./ 2.^(0:ceil(log2(top/min(lam))));      % down to and past the finest mode
    nb = numel(edges) - 1;
    bnd = nan(nb, 1);  ctr = sqrt(edges(1:nb) .* edges(2:end))';
    for i = 1:nb
        m = lam <= edges(i) & lam > edges(i+1);
        if any(m), bnd(i) = mean(d(m)); end
    end
    y = bnd / max(bnd);
    ok = ~isnan(y);  yy = y(ok);  xx = ctr(ok);
    [~, ip] = max(yy);
    peak = xx(ip);
    b = find((1:numel(yy))' > ip & yy < 0.5, 1);
    if isempty(b) || b == 1, half = NaN;
    else, half = exp(interp1(yy([b-1 b]), log(xx([b-1 b])), 0.5));
    end
    % ⭐ the robust summary: the gain-weighted centroid in LOG wavelength. Unlike a threshold
    % crossing it cannot jump, so it is the one to quote when comparing curves.
    w = max(yy, 0);
    centroid = exp(sum(w .* log(xx)) / max(sum(w), eps));
end

function N = i_normals(opts, nH)
% The surface's own vertex normals when they were passed (rheome.load.bases -> B.L.S.VertNormals),
% otherwise area-weighted from the faces.
    if ~isempty(opts.Normals)
        N = opts.Normals;
        if size(N,1) ~= nH
            error('inverse:resolution:normals', 'Normals is %s, expected %d rows.', ...
                  mat2str(size(N)), nH);
        end
    else
        P = opts.Vertices;  F = opts.Faces;
        if isempty(P) || isempty(F)
            error('inverse:resolution:normals', 'Pass Normals, or Vertices and Faces, for the MTF.');
        end
        fn = cross(P(F(:,2),:) - P(F(:,1),:), P(F(:,3),:) - P(F(:,1),:), 2);
        N = zeros(nH, 3);
        for c = 1:3
            N(:,c) = accumarray(F(:), repmat(fn(:,c), 3, 1), [nH 1]);
        end
    end
    N = N ./ max(vecnorm(N, 2, 2), eps);
end
% Author: Diellor Basha, 2026
