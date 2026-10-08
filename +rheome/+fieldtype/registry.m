function T = registry()
% FIELDTYPE.REGISTRY  Every kind of field this system can carry, as data.
%
%   T = rheome.fieldtype.registry()
%
% ⭐ A FIELD IS A TYPED PAYLOAD OVER A DOMAIN, and its type is four things: the element
% RANK it sits on, the BUNDLE (what each element carries), the number FIELD (real, complex,
% quaternion) and, for differential forms, the DEGREE. Those four decide the scalar count,
% decide whether an operator accepts it, and are what an operator's arrow is written in.
%
% ⚠ THE NAMES MIRROR nxr-compute/include/nxr/field_registry.h DELIBERATELY. The C++ side
% already catalogues these variants with the same ids and the same enum spellings; keeping
% them identical means a future binding needs no translation table, and it means this
% registry can be checked against that one. `status` says which variants this MATLAB tree
% actually produces or consumes today, so the table stays honest about the gap.
%
% ⚠ NO TIME AXIS, BY DESIGN. A field type is spatial. A [nV x nT] matrix is nT fields of
% type scalarVertex, and the time axis is a separate domain (rheome.domain.index) that composes
% over this one -- which is why the same type serves one frame and a whole recording.
%
% COLUMNS
%   id            the variant name, identical to the C++ registry where one exists
%   domain_kind   which domain kinds can carry it
%   rank          vertex | edge | face | element
%   bundle        scalar | tangent | ambient | immersion
%   field         real | complex | quaternion
%   n_form        na | zero | one | two           (differential degree, when it is a form)
%   repr          na | world | local_frame | intrinsic_complex | quaternion_interleaved
%   components    scalars per element
%   status        built (this tree makes or takes it) | planned (C++ has it, here it does not)
%   notes         where it turns up
%
% See also: rheome.fieldtype.describe, rheome.fieldtype.validate, rheome.operators.registry, rheome.domain.of
%
% Author: Diellor Basha, 2026

    r = {
    % id                    domains                 rank      bundle      field       n_form repr                      comp status   notes
    'scalarVertex',        {'complex','graph','points'}, 'vertex', 'scalar',   'real',      'zero','na',                    1, 'built',  'vorticity, divergence, curvature: [nV x nT]'
    'scalarFace',          {'complex'},            'face',   'scalar',   'real',      'na',  'na',                    1, 'built',  'per-face wavenumber magnitude, winding number'
    'oneFormEdge',         {'complex','graph'},    'edge',   'scalar',   'real',      'one', 'na',                    1, 'planned','no d0/d1 here; nxr-compute has them'
    'twoFormFace',         {'complex'},            'face',   'scalar',   'real',      'two', 'na',                    1, 'planned','needs a Hodge star'
    'tangentVertex',       {'complex'},            'vertex', 'tangent',  'complex',   'na',  'intrinsic_complex',     1, 'built',  'rheome.operators.connection_laplacian'
    'tangentFace',         {'complex'},            'face',   'tangent',  'complex',   'na',  'intrinsic_complex',     1, 'planned',''
    'tangentEdge',         {'complex'},            'edge',   'tangent',  'complex',   'na',  'intrinsic_complex',     1, 'planned',''
    'ambientVertexWorld',  {'complex','graph','points'}, 'vertex', 'ambient',  'real',      'na',  'world',                 3, 'built',  'current J [3nV x nT] interleaved; leadfield columns'
    'ambientVertexLocal',  {'complex'},            'vertex', 'ambient',  'real',      'na',  'local_frame',           3, 'planned','the frame lift lives in nxr-compute'
    'ambientFaceWorld',    {'complex'},            'face',   'ambient',  'real',      'na',  'world',                 3, 'built',  'per-face gradient / wavevector [nF x 3]'
    'ambientEdge',         {'complex'},            'edge',   'ambient',  'real',      'na',  'world',                 3, 'planned',''
    'immersionVertex',     {'complex'},            'vertex', 'immersion','quaternion','na',  'quaternion_interleaved',4, 'built',  'relative Dirac [4nV x 1], order [w x y z]'
    'immersionFace',       {'complex'},            'face',   'immersion','quaternion','na',  'quaternion_interleaved',4, 'planned','diracFace in nxr-compute'
    % --- this repo's additions: the spectral line and the sensor array as carriers ---
    'coeffScalar',         {'index'},              'element','scalar',   'real',      'na',  'na',                    1, 'built',  'Laplace-Beltrami coefficients [K x nT]'
    'coeffCurrent',        {'index'},              'element','scalar',   'real',      'na',  'na',                    1, 'built',  'relative-Dirac modal coefficients of a current: one scalar per mode'
    'sensorScalar',        {'points','complex'},   'vertex', 'scalar',   'real',      'na',  'na',                    1, 'built',  'the measured field [nCh x nT]'
    };
    T = table(string(r(:,1)), r(:,2), string(r(:,3)), string(r(:,4)), string(r(:,5)), ...
              string(r(:,6)), string(r(:,7)), cell2mat(r(:,8)), string(r(:,9)), string(r(:,10)), ...
              'VariableNames', {'id','domain_kind','rank','bundle','field','n_form','repr','components','status','notes'});
end
% Author: Diellor Basha, 2026
