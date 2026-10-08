function T = registry()
% OPERATORS.REGISTRY  Every operator as a typed arrow: what it takes, what it gives, its shape.
%
%   T = rheome.operators.registry()
%   T = rheome.operators.registry(); T(T.in == "ambientVertexWorld", :)     % what accepts a current
%
% ⭐ THIS IS docs/OPERATOR-REGISTRY.md MADE MACHINE-READABLE. That document is canonical for
% the NAMES; this table is canonical for the ARROWS. Every row says which field type the
% operator consumes, which it produces, and the shape of the matrix that does it as a
% function of the domain's element counts -- so a chain can be checked before it is run, and
% the shapes can be checked against the code that builds them (tests/tOperatorRegistry.m
% builds each mesh operator on an icosphere and compares).
%
% ⚠⚠ DIRECTION IS NOT THE SAME AS CONVENTION, and this is where the arrows stop protecting you.
% Two operators can share `in` and `out` and still not be interchangeable, because one is the
% ANALYSIS map and the other the SYNTHESIS map for the same pair of spaces. `forward_dirac` and
% `leadfield_dirac` are both coeffCurrent -> sensorScalar; they differ by a mass weighting of about
% 1/(mean vertex area) = 9e4, measured. Substituting one for the other type-checks, runs, returns
% tesla and produces the right spatial pattern 100+ dB too small. ⭐ The `notes` column is the only
% place that distinction lives, so read it before chaining a source-mapping arrow.
%
% ⚠⚠ AND A COEFFICIENT AXIS IS NOT A LENGTH AXIS. `coeffScalar` (Laplace-Beltrami) carries
% eigenvalues with a calibrated wavelength, lambda = 2*pi/sqrt(eig). `coeffCurrent` (Dirac) does NOT:
% its spectrum has no metric interpretation, so a filter designed in millimetres --
% rheome.filters.mexhat(Lambda, (mm/2pi)^2) -- is dimensionally meaningless on it and fails QUIETLY. Measured
% consequence: a detection sweep returned the same threshold for 267, 189, 140, 95 and 70 mm to three
% digits, because the filter landed in the same nearly-flat region of an unrelated spectrum every
% time. ⭐ Calibrate first by projecting each Dirac mode's current part onto the LBO and taking its
% energy-weighted mean wavenumber; then the filter is on a real length axis. The two field types do
% not encode this difference, so it cannot be caught by rheome.pipeline.then -- it is a caveat, not a guard.
%
% ⚠ RANK CHANGES ARE THE POINT. face_gradient takes a vertex scalar to a face vector,
% face_average takes it back, weak_curl takes an ambient vertex field to a vertex scalar.
% Read the `in`/`out` columns and the direction is never in doubt.
%
% ⚠ THIS TREE HAS NO EDGE OPERATORS. There is no d0, no d1 and no Hodge star here: the
% discretisation is FEM hat functions on faces, not DEC on a halfedge complex. The edge
% variants exist in nxr-compute and are marked 'planned' in rheome.fieldtype.registry so the gap
% is visible rather than surprising.
%
% COLUMNS
%   id       the operator, as called
%   fcn      the function that builds or applies it
%   in, out  field type ids (rheome.fieldtype.registry)
%   dom      which domain the operator is built on: 'mesh' | 'cross' (domain to domain)
%   shape    @(d) [rows cols] against a domain descriptor, or [] when it is not a matrix
%   apply    @(env, X) Y -- applies the operator to a field, given a bound environment
%            (rheome.pipeline.env). Empty where a pipeline cannot run it yet, and rheome.pipeline.compile
%            refuses such a step BY NAME rather than failing at run time.
%   target   which of the environment's domains the output lives on: '' = the same one,
%            else a key into env.domains ('source', 'modes', 'sensors')
%   status   built | planned
%   notes
%   direction   'forward' | 'inverse' | 'endo' (within one manifold) -- which way across manifolds
%   convention  'analysis' | 'synthesis' | '' -- WHICH of the two maps between the same pair of
%               spaces, because direction alone does not determine the operator (see below)
%   adjoint     the paired operator id where the pair is a true adjoint pair, '' otherwise
%
% ⚠ NOT EVERY OPERATOR IS A MATRIX, and the last three rows are the reason the `apply` column
% exists as a function handle rather than a kernel. rheome.flow.activation takes a modulus and
% rheome.dynamics.opticalflow_scalar solves a variational problem per frame: both are NONLINEAR in
% the data, so neither can be fused into the imaging kernel the way rheome.flow.curl is. They also
% change the frame axis -- optical flow consumes nT frames and returns nT-1 -- which a pure
% field-type arrow does not express; rheome.pipeline.run checks the shapes, not the frame count.
%
% See also: rheome.fieldtype.registry, rheome.domain.of, rheome.operators.check, docs/OPERATOR-REGISTRY.md
%
% Author: Diellor Basha, 2026

    r = {
    % id                     fcn                                   in                    out                   dom      shape                       apply                                                     target     status    notes
    'mass',                  'rheome.operators.mass',                     'scalarVertex',       'scalarVertex',       'mesh',  @(d) [d.nV d.nV],           @(e, X) e.M * X,                                          '',        'built',  'the inner product on vertex scalars'
    'laplace_beltrami',      'rheome.operators.laplace_beltrami',         'scalarVertex',       'scalarVertex',       'mesh',  @(d) [d.nV d.nV],           @(e, X) e.L * X,                                          '',        'built',  'cotangent stiffness; PSD, rows sum to zero'
    'face_gradient',         'rheome.operators.face_gradient',            'scalarVertex',       'ambientFaceWorld',   'mesh',  @(d) [d.nF d.nV],           @(e, X) i_interleave(e.fg.Gx*X, e.fg.Gy*X, e.fg.Gz*X),    '',        'built',  'Gx/Gy/Gz, one [nF x nV] block per ambient component'
    'face_average',          'rheome.operators.face_gradient (.W)',       'scalarFace',         'scalarVertex',       'mesh',  @(d) [d.nV d.nF],           @(e, X) e.fg.W * X,                                       '',        'built',  'area-weighted face to vertex, rows normalised'
    'weak_curl',             'rheome.operators.weak_differential (.Curl)','ambientVertexWorld', 'scalarVertex',       'mesh',  @(d) [d.nV 3*d.nV],         @(e, X) e.wd.Curl * X,                                    '',        'built',  'integrates by parts; no mass matrix afterwards'
    'weak_div',              'rheome.operators.weak_differential (.Div)', 'ambientVertexWorld', 'scalarVertex',       'mesh',  @(d) [d.nV 3*d.nV],         @(e, X) e.wd.Div * X,                                     '',        'built',  'as weak_curl, with grad psi instead of N x grad psi'
    'curl',                  'rheome.differential.curl',                  'ambientVertexWorld', 'scalarVertex',       'mesh',  [],                         @(e, X) rheome.differential.curl(X, e.S, e.fg),                  '',        'built',  'strong form: faces then fg.W; CCW positive'
    'divergence',            'rheome.differential.divergence',            'ambientVertexWorld', 'scalarVertex',       'mesh',  [],                         @(e, X) rheome.differential.divergence(X, e.S, e.fg),            '',        'built',  'strong form; source positive'
    'poisson',               'rheome.differential.poisson',               'scalarVertex',       'scalarVertex',       'mesh',  [],                         @(e, X) rheome.differential.poisson(X, e.S),                     '',        'built',  'K phi = M f'
    'connection_laplacian',  'rheome.operators.connection_laplacian (.A)','tangentVertex',      'tangentVertex',      'mesh',  @(d) [d.nV d.nV],           [],                                                       '',        'built',  'complex Hermitian; the builder returns a struct, .A is the operator'
    'dirac_extrinsic',       'rheome.operators.dirac_extrinsic',          'immersionVertex',    'immersionVertex',    'mesh',  @(d) [4*d.nV 4*d.nV],       [],                                                       '',        'built',  'already squared: E = D'' W_F D'
    'dirac_intrinsic_sq',    'rheome.operators.dirac_intrinsic_sq',       'immersionVertex',    'immersionVertex',    'mesh',  @(d) [4*d.nV 4*d.nV],       [],                                                       '',        'built',  'the intrinsic block'
    'dirac_frame',           'rheome.operators.dirac_frame',              'immersionVertex',    'immersionVertex',    'mesh',  @(d) [4*d.nV 4*d.nV],       [],                                                       '',        'built',  'the relative Dirac: tau blend of the two blocks'
    'phase_gradient',        'rheome.flow.phasegradient',                 'scalarVertex',       'ambientFaceWorld',   'mesh',  [],                         [],                                                       '',        'built',  'k per face, omega per vertex (two ranks out)'
    'phase_singularity',     'rheome.detect.phasesingularity',            'scalarVertex',       'scalarFace',         'mesh',  [],                         [],                                                       '',        'built',  'winding number per triangle'
    'd0',                    '(nxr-compute)',                      'scalarVertex',       'oneFormEdge',        'mesh',  @(d) [d.nE d.nV],           [],                                                       '',        'planned','exterior derivative; no edge operators in this tree'
    'hodge1',                '(nxr-compute)',                      'oneFormEdge',        'oneFormEdge',        'mesh',  @(d) [d.nE d.nE],           [],                                                       '',        'planned','needs the dual mesh'
    % ---- cross-domain: the maps between carriers ----
    'lb_forward',            'rheome.eigen.modes (Phi'' M)',              'scalarVertex',       'coeffScalar',        'cross', [],                         @(e, X) e.lbo.Phi' * (e.lbo.Mass * X),                    'modes',   'built',  'cortex to the eigenmode line'
    'lb_inverse',            'rheome.eigen.modes (Phi)',                  'coeffScalar',        'scalarVertex',       'cross', [],                         @(e, X) e.lbo.Phi * X,                                    'source',  'built',  'the line back to the cortex: this is the lift for drawing'
    'inverse_mne',           'rheome.inverse.mne',                        'sensorScalar',       'ambientVertexWorld', 'cross', [],                         @(e, X) e.K * X,                                          'source',  'built',  'ImagingKernel [3nV x nCh]; nComponents is its fibre tag'
    'inverse_dirac',         'rheome.inverse.dirac',                      'sensorScalar',       'coeffCurrent',       'cross', [],                         [],                                                       'modes',   'built',  'straight to modal coefficients, no vertex field'
    'forward_project',       'rheome.forward.project',                    'ambientVertexWorld', 'coeffCurrent',       'cross', [],                         [],                                                       'modes',   'built',  'the adjoint of reconstruct'
    'forward_reconstruct',   'rheome.forward.reconstruct',                'coeffCurrent',       'ambientVertexWorld', 'cross', [],                         [],                                                       'source',  'built',  'modal coefficients to a current field'
    'leadfield',             'rheome.forward.leadfield',                      'ambientVertexWorld', 'sensorScalar',       'cross', [],                         @(e, X) e.G * X,                                          'sensors', 'built',  'THE FORWARD DIRECTION OF SOURCE MAPPING, and it was missing from this table while every inverse was present. Gain [nCh x 3nV], T per A m.'
    'leadfield_dirac',       'rheome.forward.diracgain',         'coeffCurrent',       'sensorScalar',       'cross', [],                         [],                                                       'sensors', 'built',  'the SYNTHESIS gain: one column per Dirac eigenmode sensor pattern. NOT rheome.forward.dirac -- see that row.'
    'forward_dirac',         'rheome.forward.dirac',                      'coeffCurrent',       'sensorScalar',       'cross', [],                         [],                                                       'sensors', 'built',  'ANALYSIS convention: what rheome.inverse.dirac consumes. Carries a mass weighting, so it differs from leadfield_dirac by ~1/(vertex area) = 9e4. Using it to synthesise puts the field 100+ dB low, in correct units, with the right pattern shape.'
    'flow_field',            'rheome.flow.field',                         'sensorScalar',       'ambientVertexWorld', 'cross', [],                         [],                                                       'source',  'built',  'fused: sensors to per-vertex current'
    'flow_curl',             'rheome.flow.curl',                          'sensorScalar',       'coeffScalar',        'cross', [],                         [],                                                       'modes',   'built',  'fused: sensors to vorticity coefficients [Ks x C]'
    'flow_divergence',       'rheome.flow.divergence',                    'sensorScalar',       'coeffScalar',        'cross', [],                         [],                                                       'modes',   'built',  'fused: sensors to divergence coefficients'
    'flow_potential',        'rheome.flow.potential',                     'sensorScalar',       'coeffScalar',        'cross', [],                         [],                                                       'modes',   'built',  'Lambda^-1 times flow_divergence'
    'flow_stream',           'rheome.flow.stream',                        'sensorScalar',       'coeffScalar',        'cross', [],                         [],                                                       'modes',   'built',  'Lambda^-1 times flow_curl'
    % ---- pipeline 2: the KINEMATIC route. Not fusable -- these are nonlinear in the data. ----
    'activation',            'rheome.flow.activation',                    'ambientVertexWorld', 'scalarVertex',       'mesh',  [],                         @(e, X) rheome.flow.activation(X),                               '',        'built',  'current to a scalar activation map; NONLINEAR (a modulus), so no kernel'
    'optical_flow',          'rheome.dynamics.opticalflow_scalar',        'scalarVertex',       'ambientVertexWorld', 'mesh',  [],                         @(e, X) rheome.dynamics.opticalflow_scalar(X, e.S),              '',        'built',  'Horn-Schunck apparent velocity; needs >= 2 frames and LOSES one'
    'apparent_flow',         'rheome.flow.apparent',                      'scalarVertex',       'scalarVertex',       'mesh',  [],                         [],                                                       '',        'built',  'activation map to div v / curl v / |v|; a bundle, not a single field'
    };
    T = table(string(r(:,1)), string(r(:,2)), string(r(:,3)), string(r(:,4)), string(r(:,5)), ...
              r(:,6), r(:,7), string(r(:,8)), string(r(:,9)), string(r(:,10)), ...
              'VariableNames', {'id','fcn','in','out','dom','shape','apply','target','status','notes'});
    T = i_conventions(T);
end

function T = i_conventions(T)
% ⭐ DIRECTION AND CONVENTION ARE TWO AXES, NOT ONE, and that is the whole point of these columns.
%   direction   which way across the manifolds: 'forward' | 'inverse' | 'endo' (within one)
%   convention  which of the two maps between the SAME pair of spaces: 'analysis' | 'synthesis'
%   adjoint     the paired operator's id, where the pair really is an adjoint pair
%
% ⚠⚠ A SINGLE forward/inverse FLAG WOULD NOT CATCH THE BUG THAT MOTIVATED THIS. forward_dirac and
% leadfield_dirac have the SAME arrow -- coeffCurrent -> sensorScalar -- and the same direction, and
% are not interchangeable: measured on a reference subject they differ by a factor of 1.2e5, which is
% 1/(mean vertex area) = 9.1e4. Substituting one for the other type-checks, runs, returns tesla,
% produces the right spatial pattern and is 100+ dB wrong. `convention` is what separates them, and
% tests/tSelectionRegistry asserts that any two rows sharing (in, out) differ in it.
%
% ⚠ AND AN INVERSE IS NOT ALWAYS AN ADJOINT. inverse_mne is a REGULARISED pseudo-inverse: it has a
% direction but no adjoint partner, because G'*(GG' + lambda*C)^-1 is not the adjoint of G. Its
% `adjoint` is deliberately empty. Only the eigenbasis pairs -- where synthesis is Phi and analysis is
% Phi'*M -- are true adjoints, and the test checks the pairing is symmetric and swaps in/out.
%
% ⚠⚠ AND leadfield_dirac / forward_dirac ARE NOT AN ADJOINT PAIR EITHER, though it is tempting to
% register them as one. An adjoint SWAPS in and out; these two share the same arrow
% (coeffCurrent -> sensorScalar) and differ only by a mass weighting -- they are one map under two
% scalings, not two maps. Their `adjoint` is empty and `convention` alone separates them. The true
% adjoint of leadfield_dirac would be sensorScalar -> coeffCurrent, which is inverse_dirac's ARROW but
% not its content, since that one is regularised.
%
% ⚠ Two rows may legitimately share an arrow when they are simply different computations:
% flow_curl, flow_divergence, flow_potential and flow_stream are all sensorScalar -> coeffScalar, and
% nothing is wrong with that. The convention test therefore applies only among rows that CARRY a
% convention -- an analysis/synthesis ambiguity -- not to every shared arrow.
    n = height(T);
    dirn = repmat("endo", n, 1);  conv = repmat("", n, 1);  adj = repmat("", n, 1);
    K = { ...
      'lb_forward',         'forward', 'analysis',  'lb_inverse'
      'lb_inverse',         'inverse', 'synthesis', 'lb_forward'
      'forward_project',    'forward', 'analysis',  'forward_reconstruct'
      'forward_reconstruct','inverse', 'synthesis', 'forward_project'
      'leadfield',          'forward', 'synthesis', ''
      'leadfield_dirac',    'forward', 'synthesis', ''
      'forward_dirac',      'forward', 'analysis',  ''
      'inverse_mne',        'inverse', '',          ''
      'inverse_dirac',      'inverse', '',          ''
      'flow_field',         'inverse', '',          ''
      'flow_curl',          'inverse', '',          ''
      'flow_divergence',    'inverse', '',          ''
      'flow_potential',     'inverse', '',          ''
      'flow_stream',        'inverse', '',          ''
    };
    for i = 1:size(K,1)
        j = T.id == string(K{i,1});
        if any(j), dirn(j) = string(K{i,2});  conv(j) = string(K{i,3});  adj(j) = string(K{i,4}); end
    end
    T.direction = dirn;  T.convention = conv;  T.adjoint = adj;
end

% An ambient field is stored interleaved, rows [x1 y1 z1 x2 y2 z2 ...] -- the layout the
% weak operators and the imaging kernels already use, so a face gradient composes with them
% without a reshape at the call site.
function Y = i_interleave(X, Yc, Z)
    [n, m] = size(X);
    Y = zeros(3*n, m);
    Y(1:3:end, :) = X;  Y(2:3:end, :) = Yc;  Y(3:3:end, :) = Z;
end
% Author: Diellor Basha, 2026
