function out = critical(dbasis, aCoef, S, frames, opts)
% FLOW.CRITICAL  Classified critical points, with topological charge, at selected frames.
%
%   out = rheome.flow.critical(dbasis, aCoef, S, frames)
%   out = rheome.flow.critical(dbasis, aCoef, S, frames, opts)
%
% THE VECTOR BRANCH. Everything else in this pipeline works on SCALARS -- curl gives vortices,
% divergence gives sources and sinks -- and scalars cannot carry a topological charge. A critical
% point is a ZERO OF THE VECTOR FIELD, and its charge is the winding number of the complex tangent
% field around it. That needs J itself, so this is the one place the [3N] field is reconstructed.
%
% WHAT THE SCALAR BRANCH MISSES. Sources and sinks have div ~= 0 but curl ~ 0, so a curl detector
% never sees them. Saddles have BOTH ~ 0, so no scalar sees them at all. And a local maximum of a
% scalar has no index, so Poincare-Hopf says nothing about it. Only here is sum(index) = chi
% testable.
%
% COST. J is [3N x nFrames], so this runs on SELECTED frames, not all of them -- typically the
% high-amplitude windows found with rheome.show.activity, whose WHEN panel is free. The Dirac
% coefficients stay [P x nT] throughout; only the requested frames are ever expanded.
%
% INPUTS:
%   dbasis  Dirac eigenbasis (.Phi, .nVert, .nModes)
%   aCoef   [P x nT] Dirac mode coefficients over time (rheome.flow.synth on the Dirac joint spectrum)
%   S       surface (.Vertices, .Faces, .Hemi)
%   frames  indices into aCoef columns
%   opts    .types  passed to rheome.detect.criticalPoints (default 'all')
%           .op     a precomputed rheome.detect.operator(S); built once here if omitted
%
% OUTPUT (struct out):
%   .cp     [1 x nFrames] struct array from rheome.detect.criticalPoints, one per frame
%   .frames the frame indices used
%   .counts [nFrames x 4] vortex / source / sink / saddle per frame
%   .chi    [nFrames x nHemi] index sum per hemisphere -- the Poincare-Hopf check
%   .chiExpected  Euler characteristic per hemisphere (2 for a closed surface)
%
% ⚠ .chi is a TEST, not a result. It equals the Euler characteristic only if the mesh is closed
% and manifold; a mismatch indicates mesh trouble (non-manifold edges, boundary) rather than
% anything about the data. Report it either way.
%
% See also: rheome.detect.criticalPoints, rheome.detect.operator, rheome.forward.reconstruct, rheome.flow.observables
%
% Author: Diellor Basha, 2026

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'types') || isempty(opts.types), opts.types = 'all'; end
    if ~isfield(opts,'op')    || isempty(opts.op),    opts.op    = rheome.detect.operator(S); end
    frames = frames(:).';

    J = rheome.forward.reconstruct(aCoef(:, frames), dbasis);      % [3N x nFrames] -- only these frames
    if ~isreal(J), J = real(J); end                          % critical points of the real field

    nF = numel(frames);
    cp = cell(1, nF);  counts = zeros(nF, 4);  chi = [];
    kinds = {'vortex','source','sink','saddle'};
    for k = 1:nF
        c = rheome.detect.criticalPoints(J(:,k), S, opts.types, opts.op);
        cp{k} = c;
        for i = 1:4, counts(k,i) = sum(strcmpi(c.type, kinds{i})); end
        if isfield(c,'chi')
            if isempty(chi), chi = zeros(nF, numel(c.chi)); end
            chi(k,:) = c.chi(:).';                            %#ok<AGROW>
        end
    end

    out.cp     = [cp{:}];
    out.frames = frames;
    out.counts = counts;
    out.kinds  = kinds;
    out.chi    = chi;
    out.chiExpected = 2;                                      % closed genus-0 surface per hemisphere
end

% Author: Diellor Basha, 2026
