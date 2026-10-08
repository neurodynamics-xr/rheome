function bank = vortexbank(name, varargin)
% FLOW.VORTEXBANK  A dictionary of vortex atoms in SENSOR space, for decomposing a band.
%
%   bank = rheome.flow.vortexbank('sub01')
%   bank = rheome.flow.vortexbank(name, Locations=80, Scales=[70 140 267], Band=[8 16])
%
% Forwards a grid of rheome.flow.vortexatom spatial factors through the leadfield once, so that matching
% a record costs two [nCh x nT] products and no per-atom synthesis. Each atom keeps its rank-two
% structure: the sensor signature is pA*a(t) + pB*b(t), and rheome.flow.vortexmatch marginalises the
% temporal phase analytically rather than by searching it.
%
% ⭐⭐ WHAT THIS DICTIONARY CAN AND CANNOT SETTLE, measured at SNR 3 against real empty-room noise
% (docs/2026-09-26-feature-table-design.md section 32):
%   LOCATION   ⭐ grid-limited, not noise-limited. An off-grid vortex is matched to the nearest
%              atom with 0.0 mm of excess error over the grid's own best, at 31 mm spacing. Refine
%              the grid and the answer refines with it.
%   PRESENCE   ⭐ decisive. Best match energy is 0.90 of the data energy for a planted vortex and
%              0.0001 for empty-room noise alone -- a ratio of 9299.
%   SCALE      ⚠ nearly degenerate. Atoms at ONE location but different scales have mutual
%              coherence 0.840, against 0.211 for different locations at one scale, and a planted
%              100 mm vortex is reported as 70 mm 40% of the time and 140 mm 60%. It brackets; it
%              does not resolve.
%   ROTATION   ⚠⚠ NOT SETTLED BY MATCH ENERGY. Over 40 locations the rotating field scores a very
%              tight 0.900 (range 0.899-0.901) while a STANDING one scores 0.697 with a range of
%              0.232 to 0.885, so the per-location ratio runs 1.02 to 3.87 with a median of 1.29
%              and 28 of 40 locations below 1.5. Any focal alpha field projects heavily onto the
%              nearest atom's two-dimensional span, so high match energy means focal and located,
%              NOT rotating. ⚠ Because the spread is so wide, a single seed proves nothing either
%              way. ⭐ Use rheome.flow.vortexmatch's .icoh, which separated rotating from standing at
%              every one of the 40 locations -- though by less than the headline suggests: worst
%              rotating 0.080 against best standing 0.023.
%
% ⚠ Phase-invariant matching cannot separate a source from a vortex, and this is structural rather
%   than a choice of statistic: at winding 1 an in-plane rotation IS a temporal phase shift
%   (section 30), so the Helmholtz phase and the temporal phase are the same dial. Marginalising
%   one marginalises the other. Only the SIGN of the rotation escapes, which is why .icoh exists.
% ⚠ And the sign escapes only if the location is pinned: on a coarse grid .icoh's sign is right
%   62% of the time and wrong DETERMINISTICALLY per site, because an off-grid atom's quadrature
%   frame is rotated against the truth's (section 54). Handedness is present in the data -- exact
%   templates classify it 100% at SNR 3, with cos(B+,B-) a median -0.225 -- so what loses it is
%   this dictionary's spatial resolution, not the array.
%
% INPUTS
%   Locations  number of seed vertices by farthest-point sampling, or an explicit vector (40)
%   Scales     wavelengths in mm ([70 140 267])
%   Winding    m (1)                         Band  [8 16]
%   Duration   seconds (4)                   SampleRate  Hz (300)
%   Hemi ("L")   Bases   Gauge   Study (a rheome.load.study output)   Verbose (true)
%
% OUTPUT (struct bank)
%   .PA .PB  [nCh x N]  the forwarded quadrature pair per atom
%   .a .b    [1 x nT]   the shared temporal factor     .n0 .n1  [N x 1] atom norms
%   .vertex .scaleMM [N x 1]   .fc .band .nT .fs .chan
%
% See also: rheome.flow.vortexmatch, rheome.flow.vortexatom, rheome.flow.seedvortex, rheome.forward.leadfield
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Locations', 40);
    p.addParameter('Scales', [70 140 267], @isnumeric);
    p.addParameter('Winding', 1, @(x) isscalar(x) && x == fix(x));
    p.addParameter('Band', [8 16], @(x) numel(x)==2 && x(2)>x(1));
    p.addParameter('Duration', 4, @isscalar);
    p.addParameter('SampleRate', 300, @isscalar);
    p.addParameter('Hemi', "L");
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Gauge', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Study', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Verbose', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    B = o.Bases;  if isempty(B), B = rheome.load.bases(name); end
    Hm = B.(char(o.Hemi));  S = Hm.S;  V = S.Vertices;  nV = size(V,1);
    g = o.Gauge;
    if isempty(g), g = rheome.operators.gauge(V, double(S.Faces), Method="diffusion"); end
    st = o.Study;  if isempty(st), st = rheome.load.study(name); end
    [G, ~, chan] = rheome.forward.leadfield(st, GlobalVertices=double(Hm.gv(:))');

    % the shared temporal factor, taken from a single atom so the two agree by construction
    base = rheome.flow.vortexatom(name, WavelengthMM=o.Scales(1), Band=o.Band, Duration=o.Duration, ...
               SampleRate=o.SampleRate, Winding=o.Winding, Hemi=o.Hemi, Bases=B, Gauge=g, ...
               Check=false);

    % ⚠ farthest-point in AMBIENT space, not geodesic: a geodesic solve per candidate would cost
    %   more than the whole dictionary. The grid is a sampling, not a metric claim.
    if isscalar(o.Locations)
        nLoc = round(o.Locations);
        sel = zeros(1, nLoc);  sel(1) = round(nV/2);
        dmin = pdist2(V, V(sel(1),:));
        for i = 2:nLoc
            [~, sel(i)] = max(dmin);
            dmin = min(dmin, pdist2(V, V(sel(i),:)));
        end
    else
        sel = o.Locations(:)';  nLoc = numel(sel);
    end

    N  = nLoc*numel(o.Scales);
    PA = zeros(size(G,1), N);  PB = zeros(size(G,1), N);
    vx = zeros(N,1);  sc = zeros(N,1);  k = 0;
    for iL = 1:nLoc
        for iS = 1:numel(o.Scales)
            k = k + 1;  vx(k) = sel(iL);  sc(k) = o.Scales(iS);
            sv = rheome.flow.seedvortex(name, Vertex=sel(iL), WavelengthMM=o.Scales(iS), ...
                     Winding=o.Winding, Hemi=o.Hemi, Bases=B, Gauge=g, Check=false);
            z  = sv.z;
            nr = max(vecnorm(real(z).*g.e1 + imag(z).*g.e2, 2, 2));
            PA(:,k) = G * (reshape((  real(z).*g.e1 + imag(z).*g.e2 )', [], 1) / nr);
            PB(:,k) = G * (reshape(( -imag(z).*g.e1 + real(z).*g.e2 )', [], 1) / nr);
        end
        if o.Verbose && mod(iL, 10) == 0
            fprintf('  vortexbank: %d/%d locations\n', iL, nLoc);
        end
    end

    a = base.a;  b = base.b;
    aa = a*a';  bb = b*b';  ab = a*b';
    bank = struct('PA', PA, 'PB', PB, 'a', a, 'b', b, 'aa', aa, 'bb', bb, 'ab', ab, ...
        'n0', aa*sum(PA.^2,1)' + bb*sum(PB.^2,1)' + 2*ab*sum(PA.*PB,1)', ...
        'n1', aa*sum(PB.^2,1)' + bb*sum(PA.^2,1)' - 2*ab*sum(PA.*PB,1)', ...
        'vertex', vx, 'scaleMM', sc, 'fc', base.fc, 'band', o.Band, ...
        'nT', base.nT, 'fs', base.fs, 'chan', {chan}, 'nAtoms', N, 'locations', sel);
end

% Author: Diellor Basha, 2026
