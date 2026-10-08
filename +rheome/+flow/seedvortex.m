function out = seedvortex(name, varargin)
% FLOW.SEEDVORTEX  Seed a 3-component vortex by IMPOSING circulation in the connection frame.
%
%   out = rheome.flow.seedvortex('sub01')                       % a +1 vortex at 140 mm
%   out = rheome.flow.seedvortex(name, Winding=-1)                      % an antivortex (saddle)
%   out = rheome.flow.seedvortex(name, Phase=0)                         % a source instead
%   out = rheome.flow.seedvortex(name, Phase=pi/4)                      % a spiral
%
% ⭐⭐ WHY THIS EXISTS RATHER THAN rheome.flow.seedrotor. That one derives the field from a stream function,
% J = n x grad(psi), and the discrete co-gradient is NOT divergence-free after the face-to-vertex
% average crosses a fold -- measured div/curl = 1.26-1.34. rheome.detect.criticalPoints classifies by comparing
% |curl| against |div|, so at div ~ curl the type is a coin flip and the strongest survivors came back
% as `source` and `saddle` rather than `vortex`. Imposing the azimuthal direction directly fixes the
% TYPE CONTROL, which is the point.
%
% ⚠⚠ BUT IT DOES NOT AVOID THE CURVATURE DIVERGENCE, and I claimed it would. Measured curl/div is
% 1.06-1.23 here against 0.76 for the stream function -- a factor of ~1.5, not the order of magnitude
% implied. The reason is structural: div(A*e_azimuthal) = grad(A).e_azimuthal + A*div(e_azimuthal). The
% first term vanishes because A depends only on the geodesic distance, so grad(A) is radial and
% perpendicular to the field. The second does not: a UNIT AZIMUTHAL FIELD ON A CURVED SURFACE HAS
% DIVERGENCE, from the curvature itself. So ~1.1-1.2 is close to the floor for any vortex on this mesh,
% and no construction gets curl/div >> 1 here.
%
% THE CONSTRUCTION, in the connection Laplacian's own representation:
%   z(v) = A(d_v) * exp(i * (m*theta_v + phase))        a complex scalar per vertex
%   J(v) = Re(z)*e1(v) + Im(z)*e2(v)                    the ambient 3-vector, via the frames
% where d is the geodesic distance from the seed, theta is the angle of the outward radial direction
% measured in the gauge frame, and A vanishes at the core and decays outside it.
%
% ⭐⭐ WINDING AND TYPE ARE TWO SEPARATE DIALS, which is the thing the stream-function route could not
% give. rheome.detect.criticalPoints says it outright: "the winding is +1 for a VORTEX and a SOURCE/SINK
% alike, so only div/curl tells them apart".
%   Winding m   sets the topological charge: the direction angle advances by 2*pi*m per circuit.
%   Phase       sets what the field IS at that charge: 0 = radial (a source), pi/2 = azimuthal
%               (a vortex), in between = a spiral. m and Phase are independent.
% So a +1 source and a +1 vortex have the SAME winding and differ only by a 90 degree phase offset --
% which is why a construction that cannot control the phase cannot control the type.
%
% ⚠ THE FIELD IS PURELY TANGENTIAL by construction: it is a lift through e1/e2 and has no normal
% component. The physiologically dominant current is normal, so this is a modelling choice. It is a
% defensible one -- measured observability is -12.1 dB tangential against -12.0 dB normal -- but it is
% a choice, not a neutral default. Add a normal component explicitly if you want one.
%
% ⚠ THE SEED MUST AVOID THE GAUGE'S SINGULARITIES. The winding is m RELATIVE TO THE GAUGE, and no
% smooth frame exists on a closed surface (Poincare-Hopf forces total index chi = 2). Around a loop
% that also encloses a gauge singularity of index k, theta (the radial direction measured in the
% gauge) winds 2*pi*(1-k), so the field m*theta winds m*(1-k) in the gauge and the OBSERVED charge is
%       m + (1 - m)*k                    -- not m + k, as this note used to say.
% A +1 field reads +1 whatever k is; m = -1 at a +1 pole reads +1 (measured: a saddle seeded 2.8 mm
% from a pole of the trivial gauge came back a +1 source). `.gaugeSingularDistMM` reports how far the
% seed is from the nearest one so a bad seed is visible rather than silent.
%
% ⚠ A(d) is an ANALYTIC radial profile, (d/s)*exp(-d^2/2s^2) with s = lambda/(2*pi), not a spectral
% filter. So the scale is exact and the field is NOT band-limited to the cached basis -- which is what
% you want for a forward test (rheome.forward.simulate's Pattern route uses the vertex leadfield) and is wrong
% if you need the field to live in a mode subspace.
%
% NAME-VALUE
%   Vertex  seed (default the median-depth vertex)   WavelengthMM (140)
%   Winding m, integer (default +1)                  Phase (pi/2 = vortex; 0 = source)
%   Hemi ("L")   Gauge (a prebuilt rheome.operators.gauge)  Bases   Check (true)
%
% OUTPUT (struct out)
%   .J [3nV x 1]  .z [nV x 1] complex  .theta .dist  .vertex .winding .phase
%   .check: charge, type, distMM, nCritical, nStrong (above a 20% amplitude mask), chi, divRatio,
%           curlOverDiv, tangentialFrac, gaugeSingularDistMM
%
% ⚠ READ nStrong, NOT nCritical. Winding is undefined where |J| ~ 0 -- a singularity IS an amplitude
% zero -- so the raw count is mostly noise: 60 critical points of which 1 is strong at 70 mm, 3 at
% 140 mm and 16 at 267 mm. ⭐ A coarse pattern spans more curvature and fragments; the cleanest single
% vortex is the finest one, which is the opposite of what the instrument prefers (2*r50 = 104 mm).
%
% ⚠ rheome.detect.criticalPoints CANNOT REPORT |charge| > 1: its selection is `find(abs(q) == 1)`. So
% Winding=2 comes back as +1 and the dial is unverifiable above 1 with this detector, even if the defect
% is there.
%
% See also: rheome.flow.seedrotor, rheome.operators.gauge, rheome.operators.connection_laplacian,
%           rheome.detect.criticalPoints, rheome.forward.simulate, rheome.geom.geodesic
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Vertex', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('WavelengthMM', 140, @isscalar);
    p.addParameter('Winding', 1, @(x) isscalar(x) && x == fix(x));
    p.addParameter('Phase', pi/2, @isscalar);
    p.addParameter('Hemi', "L");
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Gauge', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Check', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    if isempty(o.Bases), B = rheome.load.bases(name); else, B = o.Bases; end
    H = B.(char(o.Hemi));  S = H.S;  nV = size(S.Vertices,1);
    V = S.Vertices;  F = double(S.Faces);
    %% the gauge: a smooth tangent frame, and where it breaks
    if isempty(o.Gauge), g = rheome.operators.gauge(V, F, Method="diffusion"); else, g = o.Gauge; end

    % ⚠⚠ THE DEFAULT SEED MUST BE CHOSEN AGAINST THE GAUGE, NOT ARBITRARILY. round(nV/2) landed 8 mm
    % from a gauge singularity on this mesh, where the observed charge is m + (1-m)*k for a gauge index k
    % and the measurement is meaningless. The default is now the vertex farthest from any of them.
    v = o.Vertex;
    if isempty(v)
        if isempty(g.singular)
            v = round(nV/2);
        else
            fcS = (V(F(g.singular,1),:) + V(F(g.singular,2),:) + V(F(g.singular,3),:))/3;
            [~, v] = max(min(pdist2(V, fcS), [], 2));
        end
    end

    %% the outward radial direction, from the gradient of the geodesic distance
    d = rheome.geom.geodesic(S, v);
    fg = rheome.operators.face_gradient(V, F);
    gx = fg.Gx*d;  gy = fg.Gy*d;  gz = fg.Gz*d;
    er = [fg.W*gx, fg.W*gy, fg.W*gz];
    er = er - sum(er.*g.normal, 2).*g.normal;                 % tangential part only
    er = er ./ max(vecnorm(er, 2, 2), eps);

    %% the angle of that direction IN THE GAUGE, and the imposed field
    theta = atan2(sum(er.*g.e2, 2), sum(er.*g.e1, 2));
    s = o.WavelengthMM*1e-3/(2*pi);
    A = (d/s) .* exp(-(d.^2)/(2*s^2));                        % 0 at the core, peak at s
    z = A .* exp(1i*(o.Winding*theta + o.Phase));
    JJ = real(z).*g.e1 + imag(z).*g.e2;
    JJ = JJ / max(vecnorm(JJ, 2, 2));
    J = reshape(JJ', [], 1);

    out = struct('J', J, 'z', z, 'theta', theta, 'dist', d, 'vertex', v, ...
                 'winding', o.Winding, 'phase', o.Phase, 'wavelengthMM', o.WavelengthMM, ...
                 'check', struct());

    %% verify
    if o.Check
        c = struct();
        if isempty(g.singular)
            c.gaugeSingularDistMM = Inf;
        else
            fc = (V(F(g.singular,1),:) + V(F(g.singular,2),:) + V(F(g.singular,3),:))/3;
            c.gaugeSingularDistMM = 1e3*min(vecnorm(fc - V(v,:), 2, 2));
        end
        dv = rheome.differential.divergence(J, S, fg);
        cv = rheome.differential.curl(J, S, fg);
        av = full(sum(H.lbo.Mass, 2));
        nrm = @(x) sqrt(sum(av.*x.^2));
        c.divRatio    = nrm(dv)/max(nrm(cv), realmin);
        c.curlOverDiv = nrm(cv)/max(nrm(dv), realmin);
        c.tangentialFrac = 1 - sum(sum(JJ.*g.normal,2).^2)/max(sum(JJ(:).^2), realmin);
        Sg = S;  Sg.nV = nV;  Sg.Hemi = {(1:nV)'};
        Sg.HemiLabel = {sprintf('Cortex %s', o.Hemi)};
        if ~isfield(Sg,'Comment') || isempty(Sg.Comment), Sg.Comment = char(name); end
        if ~isfield(Sg,'SurfaceFile'), Sg.SurfaceFile = ''; end
        cp = rheome.detect.criticalPoints(J, Sg);
        c.chi = sum(cp.chi);  c.nCritical = numel(cp.charge);
        if isempty(cp.charge)
            c.charge = NaN;  c.type = "none";  c.distMM = NaN;  c.nStrong = 0;
        else
            [~, nn] = min(pdist2(cp.pos, V), [], 2);
            a = vecnorm(JJ, 2, 2);
            strong = a(nn) > 0.20*max(a);
            c.nStrong = sum(strong);
            dd = vecnorm(cp.pos - V(v,:), 2, 2);
            [~, k] = min(dd);
            c.charge = cp.charge(k);  c.type = string(cp.type{k});  c.distMM = 1e3*dd(k);
        end
        out.check = c;
    end
end

% Author: Diellor Basha, 2026
