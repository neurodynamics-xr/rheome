function out = movingvortex(name, varargin)
% FLOW.MOVINGVORTEX  A vortex atom whose core TRANVELS, at a chosen speed along a geodesic.
%
%   out = rheome.flow.movingvortex('sub01')                   % alpha, 0.05 m/s
%   out = rheome.flow.movingvortex(name, SpeedMS=0.2)                 % faster, and broader in band
%   out = rheome.flow.movingvortex(name, SpeedMS=0, ...)              % degenerates to rheome.flow.vortexatom
%
% ⭐⭐ WHY THIS IS A SEPARATE FUNCTION. rheome.flow.vortexatom is a SEPARABLE product -- one spatial factor
% times one temporal factor -- so its core cannot move: scale, rate and place are independent by
% construction and the atom has exactly zero speed. Motion is the one thing separability forbids.
% Here the spatial factor depends on t, and the object stops being rank two.
%
% ⭐ IT IS RANK 2K, AND STILL TWO MATRIX PRODUCTS. The core is carried through K waypoints along a
% geodesic, each waypoint a vortex of the same scale and charge, blended by a time partition of unity
% W:
%       J(:,t) = JA*(W.*a)(:,t) + JB*(W.*b)(:,t),      JA, JB  [3nV x K],  W [K x nT]
% so the forward is (G*JA)*(W.*a) + (G*JB)*(W.*b): 2K leadfield columns for the whole record, and
% [3nV x nT] is still never formed.
%
% ⚠⚠ SPEED FIGHTS THE BAND, AND THE LIMIT IS TIGHT. A vertex sees the vortex only while it passes,
% for about lambda/v seconds, and that passage envelope multiplies the carrier -- so motion broadens
% the spectrum exactly as an amplitude modulation does. Staying inside the octave needs the passage to
% be no shorter than the tile's own support, v <= lambda/T_support, which at 140 mm and alpha's 1.26 s
% is about 0.11 m/s. ⭐ MEASURED, that estimate is conservative: in-band energy is 1.0000 / 0.9997 /
% 0.9984 / 0.9765 / 0.7845 / 0.4637 at 0.02 / 0.05 / 0.11 / 0.25 / 0.50 / 1.00 m/s, so the knee is
% nearer 0.5 m/s than 0.11. But the conclusion survives: at 1 m/s, a figure commonly quoted for
% cortical travelling waves, less than half the energy is still in the octave and the peak has
% dragged from 11.25 Hz down to 9.75 Hz. A fast travelling vortex is not a narrowband object.
% `.check.inbandFrac` measures it; do not assume it.
%
% ⚠⚠ TWO SPEEDS, DIFFERING BY THE FOLDING, AND THE LITERATURE'S FIGURE IS AMBIGUOUS BETWEEN THEM.
% SpeedMS is GEODESIC -- along the cortical sheet, which is what a wave in the tissue travels. The
% centroid of the field displaces through SPACE, which is what an observer outside the head would
% infer, and the two differ by however much the cortex folds in between: measured geodesic/chord
% 1.69 to 2.27 on this mesh, so a core covering 146.6 mm of surface displaces only 86.8 mm of space.
% `.check.speedGeodesicMS` and `.check.speedEuclideanMS` report both and `.check.foldFactor` their
% ratio. ⚠ Do not compare a geodesic speed here against a metres-per-second figure from a paper
% without checking which one that paper measured.
%
% ⚠ The path length is speed x the envelope support, so a fast vortex needs a long path and this mesh
% is only ~180 mm across. The path is clamped to the farthest reachable vertex and a warning fires.
%
% INPUTS
%   SpeedMS   core speed along the surface, m/s (default 0.05)
%   Vertex    the START vertex (default rheome.flow.seedvortex's gauge-aware choice)
%   Target    end vertex, overriding SpeedMS-derived distance ([] by default)
%   Waypoints K, or [] to space them at WavelengthMM/4 (default)
%   plus WavelengthMM, Winding, Phase, Band, Duration, SampleRate, Chirality, MomentNAm,
%   Hemi, Bases, Gauge -- all as rheome.flow.vortexatom
%
% OUTPUT (struct out)
%   .JA .JB [3nV x K]   .W [K x nT]   .a .b .psi   .path [1 x K] vertices
%   .pathLengthMM .speedMS .transitSec   .materialise  @() [3nV x nT]
%   .check  speedGeodesicMS, speedEuclideanMS, foldFactor, chordMM, inbandFrac, peakHz,
%           passageSec, nWaypoints, clamped
%
% See also: rheome.flow.vortexatom, rheome.flow.seedvortex, rheome.geom.geodesic, rheome.filters.travwave, rheome.dynamics.dispersion
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('SpeedMS', 0.05, @(x) isscalar(x) && x >= 0);
    p.addParameter('Vertex', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Target', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Waypoints', [], @(x) isempty(x) || (isscalar(x) && x >= 1));
    p.addParameter('WavelengthMM', 140, @isscalar);
    p.addParameter('Winding', 1, @(x) isscalar(x) && x == fix(x));
    p.addParameter('Phase', pi/2, @isscalar);
    p.addParameter('Band', [8 16], @(x) numel(x)==2 && x(2)>x(1));
    p.addParameter('Duration', 4, @isscalar);
    p.addParameter('SampleRate', 300, @isscalar);
    p.addParameter('CentreSec', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Chirality', 1, @(x) isscalar(x) && abs(x)==1);
    p.addParameter('MomentNAm', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('Hemi', "L");
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Gauge', [], @(x) isempty(x) || isstruct(x));
    p.parse(varargin{:});
    o = p.Results;

    B = o.Bases;  if isempty(B), B = rheome.load.bases(name); end
    Hm = B.(char(o.Hemi));  S = Hm.S;
    g = o.Gauge;
    if isempty(g)
        g = rheome.operators.gauge(S.Vertices, double(S.Faces), Method="diffusion");
    end

    %% 1. the stationary atom supplies the temporal factor and the start vertex
    base = rheome.flow.vortexatom(name, Vertex=o.Vertex, WavelengthMM=o.WavelengthMM, ...
                Winding=o.Winding, Phase=o.Phase, Band=o.Band, Duration=o.Duration, ...
                SampleRate=o.SampleRate, CentreSec=o.CentreSec, Chirality=o.Chirality, ...
                Hemi=o.Hemi, Bases=B, Gauge=g, Check=false);
    v0 = base.vertex;  nT = base.nT;  fs = base.fs;
    sup = base.check.supportSec;

    %% 2. the path: length = speed x the envelope's support
    dStart = rheome.geom.geodesic(S, v0);
    wantL  = o.SpeedMS * sup;                       % metres
    clamped = false;
    if isempty(o.Target)
        if wantL <= 0
            vT = v0;  L = 0;
        else
            [maxD, vFar] = max(dStart);
            if wantL > maxD
                clamped = true;
                warning('flow:movingvortex:clamped', ...
                    ['%.0f mm of travel was asked for (%.3f m/s over a %.2f s support) but the ' ...
                     'farthest vertex is %.0f mm away. Clamping; the realised speed is %.3f m/s.'], ...
                    wantL*1000, o.SpeedMS, sup, maxD*1000, maxD/sup);
                vT = vFar;  L = maxD;
            else
                [~, vT] = min(abs(dStart - wantL));  L = dStart(vT);
            end
        end
    else
        vT = o.Target;  L = dStart(vT);
    end
    dEnd = rheome.geom.geodesic(S, vT);

    % waypoints ON the geodesic: a vertex is on it when d_start + d_end is minimal
    if L == 0
        K = 1;  path = v0;  sFrac = 0;
    else
        K = o.Waypoints;
        if isempty(K), K = max(2, round(L/(o.WavelengthMM/1000/4)) + 1); end
        K = round(K);
        sFrac = linspace(0, 1, K);
        path = zeros(1, K);
        for k = 1:K
            cost = abs(dStart - sFrac(k)*L) + abs(dEnd - (1 - sFrac(k))*L);
            [~, path(k)] = min(cost);
        end
    end

    %% 3. one vortex per waypoint, sharing scale, charge and phase
    nV = size(S.Vertices, 1);
    JA = zeros(3*nV, K);  JB = zeros(3*nV, K);
    for k = 1:K
        sv = rheome.flow.seedvortex(name, Vertex=path(k), WavelengthMM=o.WavelengthMM, ...
                 Winding=o.Winding, Phase=o.Phase, Hemi=o.Hemi, Bases=B, Gauge=g, Check=false);
        zk  = sv.z;
        nrm = max(vecnorm(real(zk).*g.e1 + imag(zk).*g.e2, 2, 2));
        JA(:,k) = reshape((  real(zk).*g.e1 + imag(zk).*g.e2 )', [], 1) / nrm;
        JB(:,k) = reshape(( -imag(zk).*g.e1 + real(zk).*g.e2 )', [], 1) / nrm;
    end

    %% 4. the time partition of unity that carries the core along the path
    tt = (0:nT-1)/fs;
    [~, kPk] = max(abs(base.psi));
    tMid = tt(kPk);
    transit = L / max(o.SpeedMS, realmin);
    if L == 0 || ~isfinite(transit), transit = sup; end
    if K == 1
        W = ones(1, nT);
    else
        tk = tMid + (sFrac - 0.5)*transit;           % the core is at waypoint k at time tk(k)
        W  = zeros(K, nT);
        for k = 1:K                                   % triangular hats, summing to one
            W(k,:) = max(0, 1 - abs(tt - tk(k)) / (transit/(K-1)));
        end
        % ⚠⚠ HOLD AT THE ENDS. Before the first waypoint's time and after the last, every hat is
        %   zero, so a plain normalise leaves the field VANISHING there rather than parked. When the
        %   transit is shorter than the envelope -- which is exactly what happens once the path is
        %   clamped -- that zeroing sits INSIDE the support and broadens the spectrum all by itself,
        %   which is an artefact of the weighting and not a property of motion.
        W(1, tt < tk(1))   = 1;
        W(K, tt > tk(end)) = 1;
        cs = sum(W, 1);  cs(cs == 0) = 1;
        W  = W ./ cs;
    end

    if ~isempty(o.MomentNAm)
        amp = o.MomentNAm*1e-9 / sum(vecnorm(reshape(JA(:,1),3,[])', 2, 2));
        JA = JA*amp;  JB = JB*amp;
    end

    a = base.a;  b = base.b;
    out = struct('JA', JA, 'JB', JB, 'W', W, 'a', a, 'b', b, 'psi', base.psi, ...
        'path', path, 'sFrac', sFrac, 'pathLengthMM', L*1000, 'speedMS', L/transit, ...
        'transitSec', transit, 'fc', base.fc, 'band', o.Band, 'nT', nT, 'fs', fs, ...
        'vertex', v0, 'target', vT, 'base', base, 'check', struct());
    out.materialise = @() JA*(W.*a) + JB*(W.*b);

    %% 5. verify: the speed it actually moves at, and what motion costs in band
    Wa = W.*a;  Wb = W.*b;
    nVv = size(S.Vertices, 1);

    % per-waypoint energy map, so the centroid track costs [nV x K] rather than [3nV x nT]
    Mk = zeros(nVv, K);
    for k = 1:K, Mk(:,k) = vecnorm(reshape(JA(:,k), 3, [])', 2, 2).^2; end
    Ew  = Mk * (Wa.^2 + Wb.^2);                       % [nV x nT] energy per vertex per instant
    pw  = sum(Ew, 1);
    live = pw > 0.05*max(pw);
    cen = (Ew' * S.Vertices) ./ max(pw', realmin);    % [nT x 3] energy centroid

    % ⚠ fit over the TRANSIT, not the whole live window. Once the path is clamped the core is
    %   parked at an endpoint for most of the envelope, and a fit across all of it reports a speed
    %   several times too low -- 0.194 m/s for a 1.0 m/s geodesic request, against 0.59 expected
    %   from the fold factor.
    if K > 1
        move = live & tt >= tk(1) & tt <= tk(end);
    else
        move = false(1, nT);
    end
    if sum(move) > 2
        dtr = [0; cumsum(vecnorm(diff(cen(move,:)), 2, 2))];
        pf  = polyfit(tt(move)', dtr, 1);
        out.check.speedEuclideanMS = pf(1);
    else
        out.check.speedEuclideanMS = 0;
    end
    out.check.speedGeodesicMS = out.speedMS;
    chord = norm(S.Vertices(path(end),:) - S.Vertices(path(1),:));
    out.check.chordMM    = chord*1000;
    out.check.foldFactor = L / max(chord, realmin);
    if o.SpeedMS > 0
        out.check.passageSec = (o.WavelengthMM/1000) / o.SpeedMS;
    else
        out.check.passageSec = Inf;          % ⚠ not 1/realmin, which prints as 6.3e276
    end
    out.check.nWaypoints = K;
    out.check.clamped    = clamped;

    % the band, read on the x component of the strongest vertices at mid-path
    [~, vs] = maxk(Mk(:, max(1, round(K/2))), 200);
    Xs = zeros(numel(vs), nT);
    for k = 1:K
        A3 = reshape(JA(:,k), 3, [])';  B3 = reshape(JB(:,k), 3, [])';
        Xs = Xs + A3(vs,1)*Wa(k,:) + B3(vs,1)*Wb(k,:);
    end
    P   = abs(fft(Xs, [], 2)).^2;
    ff  = (0:nT-1)*fs/nT;  kk = 1:floor(nT/2);
    inb = ff(kk) >= o.Band(1) & ff(kk) <= o.Band(2);
    out.check.inbandFrac = median(sum(P(:,kk(inb)),2) ./ max(sum(P(:,kk),2), realmin));
    [~, ip] = max(P(:,kk), [], 2);
    out.check.peakHz = median(ff(kk(ip)));
end

% Author: Diellor Basha, 2026
