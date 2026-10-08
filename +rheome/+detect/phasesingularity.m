function out = phasesingularity(z, S)
% DETECT.PHASESINGULARITY  Topological charge of a complex SCALAR field's phase.
%
%   out = rheome.detect.phasesingularity(z, S)     % z [nV x 1] complex, one frame
%
% A point around which the phase winds through a full 2*pi is a rotational singularity. Its
% charge is that winding number, and the SIGN of the charge is the chirality -- read off the
% phase, not the amplitude.
%
%       q_face = (1/2pi) * sum over the three edges of  wrap(phi_j - phi_i)
%
% ⭐ THIS IS WHY WINDING BEATS RAW CURL. An integer cannot drift. Multiply the field by any
% positive amplitude profile -- a gain change, a different normalisation, a source that waxes
% and wanes -- and the charge and its location are untouched (tPhaseGradient asserts exactly
% this). A threshold on |curl| moves with all of them.
%
% ⚠ IT REVISES THE PREMISE IN FLOW.CRITICAL, which states that scalars cannot carry a
% topological charge. True of a REAL scalar: its zero set is a curve, not a point, and it has
% no index. A COMPLEX scalar is a map to the plane, its zeros are isolated, and they carry a
% winding number. So the analytic curl gets a charge WITHOUT reconstructing the [3N] vector
% field -- which is the cost that confines rheome.flow.critical to selected frames.
%
% NOT THE SAME OBJECT AS DETECT.VORTEX, and both are wanted. rheome.detect.vortex finds a swirl in
% the instantaneous CURRENT -- a kinematic vortex, where the flow goes round. This finds a
% point about which the TIMING winds -- where neighbouring cortex is progressively further
% through its own cycle. A travelling wave has phase singularities and no kinematic vortex;
% a standing dipole has the reverse. Their disagreement is informative, not an error.
%
% ORIENTATION IS DECIDED ONCE, GLOBALLY, AND THAT IS NOT A SHORTCUT. The winding is signed by
% the face's vertex order, which need not match the outward normal (Brainstorm's cortex meshes
% are wound opposite to theirs -- 99.2% of faces). The fix is ONE area-weighted vote over the
% whole mesh, not a per-face comparison.
%
% ⚠ FLIPPING FACES INDIVIDUALLY SILENTLY BREAKS THE TOPOLOGY. The sum over faces telescopes
% only because each interior edge is traversed once in each direction and cancels; negating one
% face of a pair destroys that. On the real cortex 346 of 40960 faces (0.84%) have an averaged
% vertex normal that disagrees with their own winding -- sharp sulcal folds, where the average
% is a bad proxy, NOT genuine inconsistency (the directed-edge test confirms the mesh is
% consistently oriented). Flipping those 346 put the net charge off zero in 100 of 401 real
% alpha frames. Measured, and the reason .total is checked rather than assumed.
%
% ⚠ THE WRAP MUST BE ANTISYMMETRIC, WHICH mod() IS NOT. wrap(x) = mod(x+pi,2pi)-pi sends BOTH
% +pi and -pi to -pi, so an edge whose phase step is exactly pi contributes the same value to
% the faces on either side instead of cancelling. That is not academic: a core sitting exactly
% antipodal to two of its triangle's vertices produces exactly that step, and the winding came
% back 0 instead of +1. atan2(sin,cos) is antisymmetric to the last bit, so interior edges
% cancel exactly and the closed-surface total is exact.
%
% A step of exactly pi is still the resolution limit -- at that point the direction the phase
% travelled is genuinely unknowable from the samples. rheome.flow.phasegradient reports .aliased for
% the same reason. Sample the field well enough that this is rare.
%
% CLOSED-SURFACE CHECK. For a complex scalar, sum(charge) = 0 on a closed surface -- every
% edge is traversed twice, in opposite senses, and cancels. (For a VECTOR field the same sum
% gives the Euler characteristic instead; that is rheome.detect.criticalPoints.) A nonzero .total on
% a closed mesh means the mesh is not closed or not consistently oriented.
%
% INPUTS:
%   z  [nV x 1] complex analytic field    S  surface (.Vertices, .Faces, .VertNormals)
%
% OUTPUT (struct out):
%   .charge     [nc x 1] nonzero winding numbers (typically +-1)
%   .pos        [nc x 3] face centroid of each
%   .chirality  [nc x 1] sign(charge)
%   .face       [nc x 1] face index
%   .total      scalar, sum over ALL faces -- 0 on a closed, consistently oriented surface
%   .perFace    [nF x 1] the full charge field
%
% See also: rheome.flow.phasegradient, rheome.detect.vortex, rheome.detect.criticalPoints
%
% Author: Diellor Basha, 2026

    if size(z, 2) ~= 1
        error('detect:phasesingularity:frame', ...
            'One frame at a time: z must be [nV x 1], got [%s].', ...
            strjoin(string(size(z)), ' x '));
    end
    if isreal(z)
        error('detect:phasesingularity:real', ...
            'z is real -- a real scalar has no winding. Apply hilbert() in time first.');
    end

    V = double(S.Vertices);  F = double(S.Faces);
    a = F(:,1);  b = F(:,2);  c = F(:,3);

    % A single global sign is only meaningful on a consistently-wound mesh. On an inconsistent
    % one the charges are individually wrong AND the total stops telescoping, so say so rather
    % than return a plausible number. (MATLAB's delaunay produces exactly such meshes.)
    dirE = [a b; b c; c a];
    [~, ~, ie] = unique(dirE, 'rows');
    if any(accumarray(ie, 1) > 1)
        error('detect:phasesingularity:orientation', ...
            ['The mesh is not consistently oriented: %d directed edges are traversed the ' ...
             'same way twice. Winding numbers need a coherent orientation -- orient the ' ...
             'faces first.'], sum(accumarray(ie,1) > 1));
    end
    phi = angle(z);
    wr  = @(x) atan2(sin(x), cos(x));      % antisymmetric -- see below

    q = ( wr(phi(b) - phi(a)) + wr(phi(c) - phi(b)) + wr(phi(a) - phi(c)) ) / (2*pi);
    q = round(q);

    % One global sign, so +1 means counter-clockwise seen from OUTSIDE. fn is unnormalised
    % (= 2*area*normal), so summing the dot products is an AREA-WEIGHTED vote: big well-formed
    % faces outvote the slivers whose averaged vertex normal is unreliable.
    fn = cross(V(b,:) - V(a,:), V(c,:) - V(a,:), 2);
    if isfield(S, 'VertNormals') && ~isempty(S.VertNormals)
        vn = double(S.VertNormals);
        on = (vn(a,:) + vn(b,:) + vn(c,:)) / 3;
    else
        on = V(a,:) + V(b,:) + V(c,:);                        % fall back to radial
    end
    if sum(sum(fn .* on, 2)) < 0, q = -q; end

    hit = find(q ~= 0);
    out.perFace   = q;
    out.total     = sum(q);
    out.charge    = q(hit);
    out.face      = hit;
    out.chirality = sign(q(hit));
    out.pos       = (V(a(hit),:) + V(b(hit),:) + V(c(hit),:)) / 3;
end

% Author: Diellor Basha, 2026
