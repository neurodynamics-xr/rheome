function e = env(opts)
% PIPELINE.ENV  The bound resources a pipeline runs against, and the domains they name.
%
%   e = rheome.pipeline.env(Surface=S, Inverse=K, Modes=basis, Sensors=chanPos)
%
% ⭐ A PIPELINE IS METADATA; AN ENVIRONMENT IS THE MATRICES. The plan says "inverse_mne then
% curl then lb_forward" and can be checked, printed and stored without any of them; running
% it needs this: the surface and its derived operators, an imaging kernel, an eigenbasis.
% Keeping them apart is what lets the same plan run on another subject.
%
% Fills in what it can derive: the face gradient and the weak operators come from the
% surface, the mass and stiffness from `rheome.operators.laplace_beltrami`, and every domain
% descriptor from the thing it describes. `.domains` is the map a cross-domain operator's
% `target` looks up.
%
% OPTIONS
%   Surface   a surface struct (Vertices/Faces) -- the source domain
%   Inverse   [3nV x nCh] imaging kernel (rheome.inverse.mne .ImagingKernel)
%   Modes     an eigenbasis struct (.Phi [nV x K], .Lambda, .Mass) -- the modes domain
%   Sensors   an [nCh x 3] position matrix or a sensor struct -- the sensors domain
%   Name      a handle for the environment (goes into the record)
%   Derive    build fg / wd / L / M from the surface (true)
%
% See also: rheome.pipeline.start, rheome.operators.registry, rheome.domain.of
%
% Author: Diellor Basha, 2026

    arguments
        opts.Surface = []
        opts.Inverse double = []
        opts.Modes   = []
        opts.Sensors = []
        opts.Name  (1,1) string = ""
        opts.Derive (1,1) logical = true
    end
    e = struct('name', char(opts.Name), 'domains', struct(), 'S', [], 'fg', [], 'wd', [], ...
               'L', [], 'M', [], 'K', [], 'lbo', [], 'built', datetime('now'));

    if ~isempty(opts.Surface)
        e.S = opts.Surface;
        e.domains.source = rheome.domain.of(e.S, Name=opts.Name + "/source");
        if opts.Derive
            V = e.S.Vertices;  F = e.S.Faces;
            e.fg = rheome.operators.face_gradient(V, F);
            e.wd = rheome.operators.weak_differential(V, F, e.fg);
            [e.L, e.M] = rheome.operators.laplace_beltrami(V, F);
        end
    end
    if ~isempty(opts.Sensors)
        s = opts.Sensors;
        if isnumeric(s), s = struct('Vertices', s, 'nV', size(s, 1)); end
        e.domains.sensors = rheome.domain.of(s, Name=opts.Name + "/sensors", Kind="points");
    end
    if ~isempty(opts.Inverse)
        e.K = opts.Inverse;
        if isfield(e.domains, 'source')
            n = rheome.domain.elements(e.domains.source, 'vertex');
            if size(e.K, 1) ~= 3*n
                error('pipeline:env:inverse', ...
                    ['The kernel has %d rows; an unconstrained kernel on this surface needs ' ...
                     '3 x %d = %d. (A constrained kernel is a different field type -- see ' ...
                     'rheome.fieldtype.registry.)'], size(e.K, 1), n, 3*n);
            end
        end
    end
    if ~isempty(opts.Modes)
        e.lbo = opts.Modes;
        K = size(e.lbo.Phi, 2);
        of = [];  if isfield(e.domains, 'source'), of = e.domains.source; end
        e.domains.modes = rheome.domain.index(K, Name=opts.Name + "/lb_modes", Of=of);
    end
end
% Author: Diellor Basha, 2026
