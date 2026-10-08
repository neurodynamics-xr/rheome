function out = chargedensity(z, S, basis, radii, varargin)
% DETECT.CHARGEDENSITY  Topological charge of a phase field, at a chosen spatial aperture.
%
%   out = rheome.detect.chargedensity(z, S, basis, radii)
%
% The per-triangle winding number (rheome.detect.phasesingularity) is the charge at the SMALLEST
% aperture a mesh admits, and at that aperture a smooth random field already carries a dense
% population of cores -- measured on real alpha, ~100 per frame, matching a phase-randomised
% surrogate in count, lifetime AND spatial distribution. A rotor has to be read at the scale
% it actually occupies, which means integrating the charge over an aperture.
%
% ⭐ NO LOOP INTEGRALS ARE NEEDED, AND THIS IS EXACT, NOT AN APPROXIMATION. By discrete Stokes
% the winding around any closed contour equals the SUM OF THE PER-FACE CHARGES it encloses --
% interior edges are traversed once each way and cancel, the same telescoping that forces the
% closed-surface total to zero. So an aperture-integrated winding is a WEIGHTED SUM of a sparse
% integer field, and a heat-kernel aperture is simply that field smoothed:
%
%       q_tau = Phi * diag(exp(-Lambda*tau)) * Phi' * c          c = charges lumped at vertices
%
% one matvec per aperture in the eigenbasis you already carry. The explicit route -- a contour
% integral per vertex per timepoint -- is ~82 million loop integrals for a single 10 s page and
% forces an ROI mask; this runs the WHOLE cortex at every aperture for a few [nV x K] products.
%
% ⚠ CHARGE IS A MEASURE, NOT A FUNCTION, so it is smoothed with Phi*E*Phi' and NOT Phi*E*Phi'*M.
% Charges are point masses; inserting the mass matrix would weight each by its vertex area and
% silently make the density depend on mesh refinement.
%
% ⚠ THE APERTURE FLOOR IS SET BY THE TRUNCATION, NOT BY THE REQUEST. With K modes the kernel
% cannot localise below the finest represented wavelength (~14 mm at K=800 on this cortex); ask
% for 5 mm and you get the 14 mm kernel back without warning. .resolved reports, per radius,
% the fraction of spectral weight retained -- below ~0.99 the aperture is truncation-limited.
%
% ⭐ THE APERTURE IS WHAT SEPARATES A ROTOR FROM THE RANDOM-FIELD BACKGROUND, and flatness
% across apertures is the sharper half of it. Measured on a 60 mm sphere at r = 26/32/40 mm:
%
%   planted unit core   q = 0.894  0.929  0.954     flatness 0.937
%   random field        max|q| =  0.96   0.93  0.78  flatness 0.528 (median)
%
% ⚠ NEITHER HALF SEPARATES ALONE AT THE TAIL, so do not threshold on one. Over 20 noise
% realisations the background's BEST vertex reached |q| = 0.932 against the core's 0.929, and
% its 99th-percentile flatness reached 0.965 against the core's 0.937. The JOINT criterion is
% what works -- measured per-vertex false rate, core accepted in every case:
%
%     |q| >= 0.8  &  flatness >= 0.85   ->  0.97%        |q| >= 0.7 & flat >= 0.85 -> 1.93%
%     |q| >= 0.8  &  flatness >= 0.90   ->  0.57%   <-- operating point
%     |q| >= 0.8  &  flatness >= 0.95   ->  0.06%   but the real core is REJECTED too
%
% 0.57% of 20484 cortical vertices is still ~117 candidates per frame; temporal persistence is
% what is expected to take that the rest of the way, not a stricter single-frame threshold.
%
% ⚠ THAT SEPARATION DEPENDS ENTIRELY ON Phi BEING MASS-ORTHONORMAL, and it fails silently
% otherwise. Phi*E*Phi' is the heat kernel only when Phi'*M*Phi = I; normalising columns
% against lumped row-sums while M is Galerkin leaves cross terms, and two cores placed
% symmetrically then read 0.29 and 5.10 instead of +-0.9. Check the basis before trusting a
% number from this function.
%
% INPUTS:
%   z      [nV x nT] complex analytic field       S  surface (.Vertices, .Faces, .VertNormals)
%   basis  struct with .Phi [nV x K], .Lambda [K x 1]   (rheome.flowpage.prepare bundle supplies both)
%   radii  [1 x nR] aperture radii in MM; heat time tau = (r/1000)^2/4
%   'Threshold'  |q| past which an aperture is counted as charged. Default 0.5.
%
% OUTPUT (struct out):
%   .q         [nV x nR x nT]  aperture-enclosed charge (a unit core reads ~1 at every radius)
%   .qmed      [nV x nT]       median over radii
%   .flatness  [nV x nT]       min|q|/max|q| over radii; ~1 = scale-invariant = a real core
%   .stable    [nV x nT]       signed count of radii agreeing past Threshold
%   .perFace   [nF x nT]       the raw integer winding field
%   .radii .tau .resolved .nCore
%
% See also: rheome.detect.phasesingularity, rheome.flow.phasegradient
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Threshold', 0.5);
    p.parse(varargin{:});
    o = p.Results;

    V = double(S.Vertices);  F = double(S.Faces);
    a = F(:,1);  b = F(:,2);  c3 = F(:,3);
    nV = size(V,1);  nT = size(z,2);  nR = numel(radii);

    dirE = [a b; b c3; c3 a];
    [~,~,ie] = unique(dirE, 'rows');
    if any(accumarray(ie,1) > 1)
        error('detect:chargedensity:orientation', ...
            'The mesh is not consistently oriented; winding numbers need a coherent one.');
    end

    phi = angle(z);
    wr  = @(x) atan2(sin(x), cos(x));                 % antisymmetric: mod() is not
    q   = round(( wr(phi(b,:) - phi(a,:)) + wr(phi(c3,:) - phi(b,:)) ...
                + wr(phi(a,:) - phi(c3,:)) ) / (2*pi));

    fn = cross(V(b,:) - V(a,:), V(c3,:) - V(a,:), 2);
    if isfield(S,'VertNormals') && ~isempty(S.VertNormals)
        vn = double(S.VertNormals);  on = (vn(a,:) + vn(b,:) + vn(c3,:))/3;
    else
        on = V(a,:) + V(b,:) + V(c3,:);
    end
    if sum(sum(fn .* on, 2)) < 0, q = -q; end          % ONE global vote -- never per face

    % Lump each face's charge onto its three vertices: a [nV x nF] incidence operator,
    % so every frame is handled by one sparse product rather than a loop over time.
    nF = size(F,1);
    P  = sparse(F(:), repmat((1:nF).', 3, 1), 1/3, nV, nF);
    Cv = P * q;                                        % [nV x nT] charge as a measure

    Phi = double(basis.Phi);  Lam = double(basis.Lambda(:));
    A   = Phi.' * Cv;                                  % [K x nT]
    P2  = Phi.^2;
    out.q = zeros(nV, nR, nT);
    out.tau = zeros(1, nR);  out.resolved = zeros(1, nR);
    for r = 1:nR
        tau = (radii(r)/1000)^2 / 4;
        e   = exp(-Lam * tau);
        out.tau(r) = tau;
        out.resolved(r) = 1 - exp(-max(Lam) * tau);   % 1 = aperture fits in the basis
        % Phi*E*Phi'*Cv is a DENSITY (charge per m^2) and would read ~1/(4*pi*tau) -- thousands
        % -- at a unit core, growing as the aperture shrinks. Dividing by the kernel's own
        % centre value K_tau(v,v) = sum_k e^{-lambda_k tau} phi_k(v)^2 turns it into a soft
        % APERTURE: a unit charge sitting at v reads exactly 1, and one at distance d reads
        % ~exp(-d^2/4tau). That is the quantity a +-2pi threshold is meant to apply to.
        kvv = P2 * e;                                  % [nV x 1] kernel at its own centre
        out.q(:,r,:) = reshape((Phi * (e .* A)) ./ max(kvv, eps), nV, 1, nT);
    end
    aq = abs(out.q);
    out.qmed     = squeeze(median(out.q, 2));
    out.flatness = squeeze(min(aq, [], 2) ./ max(max(aq, [], 2), eps));
    sg = sign(out.q) .* (aq >= o.Threshold);
    out.stable   = squeeze(sum(sg, 2));
    out.perFace = q;
    out.radii   = radii;
    out.nCore   = sum(q ~= 0, 1);
end

% Author: Diellor Basha, 2026
