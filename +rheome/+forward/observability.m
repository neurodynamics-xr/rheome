function o = observability(hm, basis, opts)
% FORWARD.OBSERVABILITY  How well the sensor array sees each spatial scale.
%
%   o = rheome.forward.observability(hm, basis)
%   o = rheome.forward.observability(hm, basis, opts)
%
% The instrument is a SPATIAL FILTER over the cortex, and this measures its gain on the same
% lambda axis the analysis filterbank lives on. For each mode it reports the whitened field a
% unit-energy current density in that mode produces:
%
%     obs(j) = || C^(-1/2) * Gc * (area .* phi_j) ||    Gc = orientation-constrained gain
%
% Whitened, so it is in units of noise standard deviations: obs(j) = 1 means that mode radiates
% exactly as much field as the noise. Plotted against each mode's spatial wavelength
% 2*pi/sqrt(lambda) it is the array's gain versus spatial scale, for this anatomy and helmet.
%
% ⚠ IT IS NOT AN "APERTURE IN mm", AND THE CURVE'S SHAPE DEPENDS ON THE BASIS YOU HAND IT.
% Read the headline below before interpreting a roll-off point as a resolution limit -- on the
% Laplace-Beltrami axis there is no clean roll-off to read, and that says more about the axis
% than about the array.
%
% ⭐ WHY THIS IS THE NUMBER THAT BOUNDS EVERY SPATIAL CLAIM. An analysis aperture applied to
% reconstructed data does not measure psi_m(x); it measures psi_m . K . Gain . x. The
% instrument's aperture is already inside and no normaliser removes it (see rheome.flow.crosstalk).
% So a spatial scale is only claimable where BOTH floors clear it: the basis floor
% (frame.SigmaFloor, set by K) and this one. For MEG the instrument floor is usually the
% binding one, and it is the one nothing in the pipeline was checking.
%
% ⚠ THE DC MODE IS NOT THE BEST-SEEN MODE. A spatially uniform current density over a closed
% cortex is very nearly a closed field and produces almost no external signal, so obs RISES
% from lambda = 0 to a peak and only then rolls off. .peakIdx marks it. Any "where does the
% curve cross -X dB" search must start from that peak and walk toward FINE wavelengths;
% starting at mode 1 finds the coarse-side rise and reports nonsense.
%
% ⭐⭐ THE HEADLINE, AND IT IS ABOUT THE BASIS, NOT THE INSTRUMENT. Read on the LAPLACE-
% BELTRAMI axis this operator looks nearly scale-blind: measured on a reference subject (270 MEG
% ch, K = 1000, whitened by the study noise covariance) observability falls only 5.6 dB from
% 250 mm to 25 mm, +1.59 dB/octave, and spatial scale explains just 18% of the variance. The
% planar argument predicts -79 dB over that range (see below), so the discrepancy is enormous.
%
% ⚠ DO NOT CONCLUDE "MEG IS SCALE-BLIND". THAT READING WAS WRONG. The LBO is INTRINSIC -- it
% is built from the metric alone and is provably blind to how the surface sits in space, so
% two isometric surfaces share its spectrum while having completely different leadfields.
% Observability is extrinsic and VECTORIAL. Measured over the same 90-293 mm window:
%
%   basis                       modes   dB/octave    R^2
%   Dirac (quaternion)            800     +10.85    0.897
%   LB unconstrained (x,y,z)     3000      +6.87    0.556
%   LB normal-constrained        1000      +4.77    0.519
%
% On a basis built for vector fields WITH the extrinsic geometry in it, scale explains 90% of
% observability. The flat LB curve measures the mismatch between a scalar intrinsic axis and
% a vector extrinsic quantity -- not a property of the array. Report the Dirac axis as primary.
%
% ⭐ WHAT FLATTENS THE LB CURVE IS FOLDING, AND THE CONTROL MATTERS. Replacing the true
% normals with a smoothed normal field restores the steep scale dependence:
%
%   orientation field      dB/octave    R^2   roll-off   mean|cos| to N   median obs
%   true normals               +1.59   0.177     5.6 dB           1.000        --
%   smoothed sigma=10 mm       +9.19   0.804    32.2 dB           0.642     -4.3 dB
%   smoothed sigma=20 mm      +11.35   0.957    39.8 dB           0.518    -13.9 dB
%   head-radial                +8.85   0.897    31.0 dB           0.482    -21.1 dB
%
% ⚠ USE SMOOTHED NORMALS FOR THAT CONTROL, NEVER HEAD-RADIAL. In a spherical conductor a
% head-RADIAL dipole is magnetically SILENT, so setting orientations radial does not "remove
% the folding" -- it also picks the direction MEG sees worst, and costs 21 dB of signal on top.
% Smoothing at sigma = 10 mm costs only 4.3 dB yet already recovers R^2 = 0.80, which is what
% isolates the flipping of the normal as the cause.
%
% ⚠ AND TANGENTIAL SOURCES ARE NOT INVISIBLE. Median observability by orientation: normal
% -12.0 dB, tangent-1 -12.1 dB, tangent-2 -15.4 dB; median tangential/normal POWER ratio 0.74.
% A vertex here is ~8 mm^2 of a smoothed mesh, so its moment is a superposition whose net
% direction need not follow the mesh normal -- which is why unconstrained and loose-orientation
% inverses are standard. Do not treat the normal constraint as physics; it is a reduction, and
% comparing a constrained basis against an unconstrained one rewards the latter for having 3x
% the degrees of freedom rather than for any geometric merit.
%
% ⚠ THREE CONSEQUENCES:
%   1. There is no single instrument aperture in mm to report beside a member's sigma. A
%      scale-based floor is the wrong shape of answer on the LB axis.
%   2. Do NOT band-limit by truncating the eigenbasis. On the LB axis the well-observed
%      subspace is spread across all 1000 modes rather than concentrated in the first 270, so
%      a mode-index cut discards well-seen fine modes and keeps poorly-seen coarse ones.
%      Project onto the SVD subspace of this operator instead. The Dirac basis concentrates
%      it better (participation ratio 448/800 against 1283/3000), though that comparison is
%      sensitive to how it is normalised by basis size and should not be leaned on hard.
%   3. What DOES bound the measurement is RANK -- 94 / 183 / 255 singular values above
%      1e-2 / 1e-3 / 1e-4 of the largest, out of 270 channels, at every scale at once. No
%      change of basis moves that number.
%
% ⚠ SOURCE CONVENTION, AND IT IS NOT A FREE CHOICE. Brainstorm's Gain columns map dipole
% MOMENTS, while an eigenmode is a current DENSITY, so the moment at a vertex is its density
% times its area -- the lumped mass diagonal. .obs uses that; .obsRaw omits it. MEASURED on
% a reference subject the two differ by a factor that itself varies 4.8x across modes, so this is
% NOT a constant rescaling that leaves the shape alone: the Brainstorm cortex is far from a
% uniform mesh and fine modes concentrate where the triangles are small. Use .obs; .obsRaw is
% kept only so the size of the discrepancy stays visible.
%
% ⚠ CHANNELS MUST MATCH THE COVARIANCE. The whitener is built from the noise covariance
% restricted to the SAME rows as the gain. A Brainstorm study carries reference, EOG and ECG
% channels whose covariance block is degenerate; including them makes C singular and the
% whitener explodes. Default is the rows opts.Channels names, and the rank of the retained
% block is reported in .rankC so a silent rank deficiency is visible.
%
% INPUTS:
%   hm     rheome.io.read.headmodel struct: .Gain [nCh x 3*nV], .GridOrient [nV x 3], .nCh, .nV
%   basis  struct with .Phi [nV x K], .Lambda [K x 1], .Mass [nV x nV]
%   opts   .Channels  row indices of hm.Gain to keep (default: all)
%          .NoiseCov  [nCh x nCh] over the SAME rows, or [] for identity whitening
%          .Reg       loading as a fraction of mean(diag(C)) (default 0.05)
%
% OUTPUT (struct o):
%   .obs         [K x 1] observability, noise sd per unit-energy mode
%   .obsRaw      [K x 1] same without the mass weighting (convention check)
%   .obsDb       [K x 1] 20*log10(obs / max(obs))
%   .wavelength  [K x 1] mm, 2*pi/sqrt(lambda)
%   .sv          [min(nCh,K) x 1] singular values of the whitened Gc*M*Phi
%   .effRank     struct with fields r1e2, r1e3, r1e4 -- modes above those relative sv cuts
%   .planarD     fitted d (m) in exp(-k*d);  .planarR2  its weighted R^2
%   .rankC .nCh .K
%
% See also: rheome.flow.crosstalk, rheome.flow.sensitivity, rheome.filters.frame, demos.inverse_spectrum
%
% Author: Diellor Basha, 2026

    if nargin < 3 || isempty(opts), opts = struct(); end
    ch  = i_opt(opts, 'Channels', []);
    C   = i_opt(opts, 'NoiseCov', []);
    reg = i_opt(opts, 'Reg', 0.05);

    G = double(hm.Gain);  nV = hm.nV;
    if isempty(ch), ch = 1:size(G,1); end
    ch = ch(:).';
    G  = G(ch, :);
    nCh = size(G,1);

    Phi = double(basis.Phi);  lam = double(basis.Lambda(:));  Mass = basis.Mass;
    K = size(Phi,2);
    if size(Phi,1) ~= nV
        error('forward:observability:size', ...
            'basis has %d vertices but the gain has %d.', size(Phi,1), nV);
    end

    % --- orientation constraint: 3 columns per vertex -> 1, along the cortical normal
    O = double(hm.GridOrient);
    if isempty(O)
        error('forward:observability:orient', ...
            'hm.GridOrient is empty; an unconstrained gain cannot be reduced to a scalar basis.');
    end
    rows = (1:3*nV).';
    cols = repelem((1:nV).', 3, 1);
    Osp  = sparse(rows, cols, reshape(O.', [], 1), 3*nV, nV);
    Gc   = G * Osp;                                        % [nCh x nV]
    clear G

    % --- whitener from the noise covariance over the SAME rows
    if isempty(C)
        W = speye(nCh);  rankC = nCh;
    else
        C = double(C);
        if ~isequal(size(C), [nCh nCh])
            error('forward:observability:cov', ...
                'NoiseCov is %s but %d channels were selected.', mat2str(size(C)), nCh);
        end
        C = (C + C.')/2;
        C = C + reg * mean(diag(C)) * eye(nCh);            % loading, so the inverse is finite
        [U,S] = eig(C, 'vector');
        rankC = sum(S > max(S) * 1e-12);
        S = max(S, max(S)*1e-12);
        W = diag(1./sqrt(S)) * U.';                        % [nCh x nCh]
    end

    area = full(sum(Mass, 2));                             % lumped vertex areas
    A    = W * (Gc * (area .* Phi));                       % [nCh x K]  moments = density x area
    Araw = W * (Gc * Phi);

    o.obs    = sqrt(sum(A.^2, 1)).';
    o.obsRaw = sqrt(sum(Araw.^2, 1)).';
    o.obsDb  = 20*log10(o.obs / max(o.obs));
    o.wavelength = 2*pi ./ sqrt(max(lam, eps)) * 1000;      % mm (vertices in m)
    o.sv     = svd(A, 'econ');
    o.effRank = struct('r1e2', sum(o.sv > o.sv(1)*1e-2), ...
                       'r1e3', sum(o.sv > o.sv(1)*1e-3), ...
                       'r1e4', sum(o.sv > o.sv(1)*1e-4));
    o.rankC = rankC;  o.nCh = nCh;  o.K = K;
    [~, o.peakIdx] = max(o.obs);                           % NOT mode 1 -- see below

    % --- fit exp(-k*d) to the measured roll-off, weighted by observability
    k  = sqrt(max(lam, eps));
    ok = isfinite(o.obs) & o.obs > 0 & k > 0;
    y  = log(o.obs(ok));  x = k(ok);  wt = o.obs(ok)/max(o.obs);
    Aw = [ones(sum(ok),1), -x] .* sqrt(wt);
    b  = Aw \ (y .* sqrt(wt));
    o.planarD = b(2);
    yh = b(1) - b(2)*x;
    o.planarR2 = 1 - sum(wt.*(y-yh).^2) / max(sum(wt.*(y-mean(y)).^2), realmin);
end

function v = i_opt(s, f, d)
    if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

% Author: Diellor Basha, 2026
