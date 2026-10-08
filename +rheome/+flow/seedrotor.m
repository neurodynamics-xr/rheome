function out = seedrotor(name, varargin)
% FLOW.SEEDROTOR  Seed a 3-component VORTEX current field of known scale at a cortical vertex.
%
%   out = rheome.flow.seedrotor('sub01')
%   out = rheome.flow.seedrotor(name, Vertex=v, WavelengthMM=140, Kernel="mexhat")
%
% ⭐⭐ THE CONSTRUCTION, AND THE ONE STEP THAT IS EASY TO GET BACKWARDS. Seed a scalar STREAM FUNCTION
% with a band-pass graph wavelet on the Laplace-Beltrami spectrum, then rotate its gradient:
%
%       psi = Phi * (g(lambda) .* (Phi' * M * delta_v))          a blob of chosen spatial scale
%       J   = n x grad(psi)                                      the vortex
%
% ⚠⚠ NOT grad(psi). A GRADIENT IS CURL-FREE BY CONSTRUCTION -- it is literally the irrotational part
% of the Helmholtz decomposition, which rheome.differential.helmholtz writes out as
%     J = grad(Phi) + N x grad(Psi) + (J.n) n + harmonic
%         irrotational  solenoidal     normal    residual
% Taking grad of a blob gives a SOURCE/SINK pair, not a rotor. The rotated (co-)gradient
% N x grad(Psi) is the solenoidal part, and Psi is the stream function whose extrema ARE the
% vortices. Same scalar, same kernel; the cross product with the normal is the whole difference.
%
% ⭐ HOW IT GIVES A VORTEX, WHICH IS NOT ABOUT THE WAVELET. J = n x grad(psi) is perpendicular to
% grad(psi) by construction, hence TANGENT TO THE CONTOURS OF psi -- the level sets are the
% streamlines. The contours around a local maximum are closed loops, so J circulates. The wavelet's
% only job is to put a clean, scale-controlled EXTREMUM at the seed; the vortex is in the critical
% point of psi, not in the kernel. That is why `heat`, a low-pass with no oscillation at all, also
% gives winding +1 -- it too has a maximum at the seed.
% ⭐ And it is why the winding is NOT a free parameter: the index of a co-gradient equals the index of
% psi's critical point. A maximum or minimum gives +1 (a centre); a saddle would give -1.
%
% ⭐ MEASURED STRUCTURE at 140 mm: psi is 1.000 at the seed with its ring minimum of -0.158 at 65 mm
% (the mexhat's negative surround), and |J| peaks at 20-30 mm where |grad psi| is largest -- the
% circulation lives on the annulus between the core and the ring, not at the centre.
%
% ⭐ WHAT THE CONSTRUCTION GUARANTEES, measured rather than assumed (see .check):
%   winding = +1    at the seed, 1-4 mm from it, at every scale and kernel tried.
% ⚠⚠ BUT IT IS NOT AN ISOLATED VORTEX, AND `.check` REPORTS ONLY THE NEAREST CRITICAL POINT, WHICH
% FLATTERED IT. Measured at 140 mm: rheome.detect.criticalPoints finds 60 critical points, charges +-1, and
% the STRONGEST is a saddle 8 mm from the seed at strength 174 against the central vortex's 93. Many of
% the 60 are low-amplitude noise -- a singularity IS an amplitude zero, so winding is undefined where
% |J| is small -- but that 8 mm saddle is not. ⭐ Read `.nCritical` alongside `.charge`, and do not
% describe this as a clean single rotor.
% ⚠ And J is tangent to psi's contours EXACTLY, while those contours are not circles on a folded
% cortex: measured over a 20-120 mm shell, |J.azimuthal|/|J| = 0.795 against |J.radial|/|J| = 0.596 in
% GEODESIC coordinates. So it is a spiral when read radially from the seed, not a concentric
% circulation. The winding is unaffected; the geometry is not what the flat-sheet picture suggests.
%   irrotational    < 1% of the energy: the field really is (almost) purely non-gradient.
% ⚠⚠ AND THE TWO I FIRST CLAIMED AND HAD TO WITHDRAW:
%   div(J) = 0      FALSE on a curved discrete surface. The AMBIENT divergence measured 1.5-1.7 times
%                   the curl before the tangential re-projection below, because the face-to-vertex
%                   average mixes faces with different normals and divergence carries curvature
%                   (div_s(f*n) = 2H*f). A co-gradient is divergence-free INTRINSICALLY and on a flat
%                   sheet, not after averaging across a fold.
%   curl = -lap(psi)  ONLY ROUGHLY, and it degrades with scale: corr 0.84 at 70 mm, 0.67 at 140 mm,
%                   0.43 at 267 mm, because a coarse pattern spans more curvature variation. Read the
%                   vorticity from rheome.differential.curl, not from -lap(psi).
%
% ⚠⚠ AND THE HELMHOLTZ BUDGET DOES NOT CLOSE, which is worth knowing before using any of these parts
% as ground truth. Measured at 267 / 140 / 70 mm, with the tangential projection on:
%     solenoidal   0.727  0.753  0.779
%     irrotational 0.008  0.006  0.006
%     normal       0.000  0.000  0.000      <- what Tangential=true buys, and it is exact
%     harmonic     0.034  0.029  0.026
%     SUM          0.770  0.789  0.810      <- NOT 1
% rheome.differential.helmholtz calls its parts orthogonal, and discretely they are orthogonal only to about
% 20% on this mesh at these scales: the energies do not add because the parts overlap with opposing
% sign. ⚠ And the surface is genus 0 (chi = 2, measured), so the harmonic space has dimension 2g = 0
% and the 3% reported as harmonic is pure discretisation rather than a real component.
%
% ⭐⭐ SO USE THIS FOR TOPOLOGY, NOT FOR MAGNITUDE. The winding is exact and is what
% rheome.detect.criticalPoints keys on, so a seeded rotor is a sound ground truth for DETECTION and for
% locating a vortex. It is not clean enough to calibrate a vorticity amplitude against.
%
% ⚠ WINDING IS NOT A FREE PARAMETER HERE. It is set by the CRITICAL POINT TYPE of psi: a maximum gives
% +1 (a centre), a saddle gives -1. A band-pass wavelet on a delta gives a maximum, so this seeds +1
% only. A -1 saddle needs a psi with a saddle -- two crossed blobs, not one -- and is not built here.
%
% ⚠⚠ THE QUATERNION ROUTE CANNOT DO THIS, and it is worth saying why rather than leaving it open.
% rheome.filters.steer(W, q0) steers a Dirac quaternion field by RIGHT multiplication with a SINGLE unit
% quaternion, which is right-H-linearity -- the very symmetry that makes the Dirac spectrum 4-fold
% degenerate. It therefore applies ONE rigid rotation to every vector in the field. A vortex needs the
% direction to rotate WITH POSITION, so no choice of q0 produces one. Steering sets a uniform dipole
% orientation (what rheome.filters.impulse's `dir` argument does); it does not make a rotor.
% ⚠ And the Dirac spectrum carries no length axis, so a wavelength-valued kernel cannot be applied to
% it without calibrating a wavenumber per mode first. The LBO route has the length axis already.
%
% ⚠ A rotor is also the pattern the array sees WORST: measured, at best about 1/3 of a rotor is
% measurable and it saturates by 20-30 mm, and location beats size at fine scales (0.187 on a gyral
% crown against 0.058 on a sulcal wall). Seeding one and forwarding it is the right experiment; expect
% the answer to be unflattering.
%
% NAME-VALUE
%   Vertex        seed vertex (local to the hemisphere); default the median-depth vertex
%   WavelengthMM  the kernel's spatial scale, default 140
%   Kernel        "mexhat" (default) | "diffgauss" | "heat"   -- ⚠ heat is low-pass, so its
%                 "vortex" is the whole hemisphere; it is allowed for comparison, not for use
%   Hemi          "L" (default)
%   Check         run the verification (default true)
%
% OUTPUT (struct out)
%   .J [3nV x 1] the vortex field   .Psi [nV x 1] the stream function   .vertex .wavelengthMM
%   .check  struct: divRatio, curlCorr, charge, type, distMM, nCritical, solenoidalFrac,
%           irrotationalFrac, normalFrac, harmFrac
%
% See also: rheome.differential.helmholtz, rheome.differential.curl, rheome.detect.criticalPoints, rheome.filters.mexhat,
%           rheome.operators.face_gradient, rheome.forward.simulate, rheome.flow.stream
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Vertex', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('WavelengthMM', 140, @isscalar);
    p.addParameter('Kernel', "mexhat");
    p.addParameter('Hemi', "L");
    p.addParameter('Check', true, @islogical);
    p.addParameter('Tangential', true, @islogical);
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.parse(varargin{:});
    o = p.Results;

    if isempty(o.Bases), B = rheome.load.bases(name); else, B = o.Bases; end
    H = B.(char(o.Hemi));  S = H.S;  lbo = H.lbo;  nV = size(S.Vertices,1);
    v = o.Vertex;  if isempty(v), v = round(nV/2); end

    %% 1. the stream function: a band-pass wavelet on the LBO, localised at the vertex
    t = (o.WavelengthMM*1e-3/(2*pi))^2;
    switch lower(string(o.Kernel))
        case "mexhat",    g = rheome.filters.mexhat(lbo.Lambda, t);
        case "diffgauss", g = rheome.filters.diffgauss(lbo.Lambda, t/2, 2*t);
        case "heat",      g = rheome.filters.heat(lbo.Lambda, t);
        otherwise, error('flow:seedrotor:kernel', 'Unknown kernel "%s".', o.Kernel);
    end
    dl = zeros(nV,1);  dl(v) = 1;
    Psi = lbo.Phi * (g(:) .* (lbo.Phi' * (lbo.Mass * dl)));
    Psi = Psi / max(abs(Psi));

    %% 2. the vortex: the ROTATED gradient, exactly as rheome.differential.helmholtz composes Vsol
    fg = rheome.operators.face_gradient(S.Vertices, double(S.Faces));
    Nf = fg.FaceNormal;
    gx = fg.Gx*Psi;  gy = fg.Gy*Psi;  gz = fg.Gz*Psi;
    Jx = fg.W*(Nf(:,2).*gz - Nf(:,3).*gy);
    Jy = fg.W*(Nf(:,3).*gx - Nf(:,1).*gz);
    Jz = fg.W*(Nf(:,1).*gy - Nf(:,2).*gx);
    % ⚠⚠ THE FACE-TO-VERTEX AVERAGE TIPS THE FIELD OUT OF THE TANGENT PLANE. Each face's
    % N_f x grad(psi) lies in ITS OWN face, and fg.W averages neighbours with different normals, so the
    % vertex field acquires a normal component -- and divergence carries curvature
    % (div_s(f*n) = 2H*f), so that component makes divergence out of nothing. Re-projecting onto the
    % vertex tangent plane is what restores the construction's guarantee.
    JJ = [Jx Jy Jz];
    if o.Tangential && ~isempty(S.VertNormals)
        Nv = S.VertNormals;  Nv = Nv ./ max(vecnorm(Nv,2,2), eps);
        JJ = JJ - sum(JJ.*Nv, 2).*Nv;
    end
    J  = reshape(JJ', [], 1);
    J  = J / max(vecnorm(JJ, 2, 2));

    out = struct('J', J, 'Psi', Psi, 'vertex', v, 'wavelengthMM', o.WavelengthMM, ...
                 'kernel', char(lower(string(o.Kernel))), 'check', struct());

    %% 3. verify the three guarantees
    if o.Check
        dv = rheome.differential.divergence(J, S, fg);
        cv = rheome.differential.curl(J, S, fg);
        lap = lbo.Phi * (lbo.Lambda(:) .* (lbo.Phi' * (lbo.Mass * Psi)));   % lap(psi), weak form
        c = struct();
        c.divRatio = norm(dv)/max(norm(cv), realmin);          % ⭐ must be << 1
        keep = abs(cv) > 0.05*max(abs(cv));                    % where the vorticity actually is
        c.curlCorr = corr(cv(keep), -lap(keep));               % ⭐ curl = -lap(psi)
        Sg = S;  Sg.nV = nV;  Sg.Hemi = {(1:nV)'};
        Sg.HemiLabel = {sprintf('Cortex %s', o.Hemi)};
        if ~isfield(Sg,'Comment') || isempty(Sg.Comment), Sg.Comment = char(name); end
        if ~isfield(Sg,'SurfaceFile'), Sg.SurfaceFile = ''; end
        cp = rheome.detect.criticalPoints(J, Sg);
        % the strongest critical point should be a vortex, at the seed
        if ~isempty(cp.charge)
            d = vecnorm(cp.pos - S.Vertices(v,:), 2, 2);
            [~, near] = min(d);
            c.charge   = cp.charge(near);
            c.type     = string(cp.type{near});
            c.distMM   = 1e3*d(near);
            c.nCritical = numel(cp.charge);
        else
            c.charge = NaN;  c.type = "none";  c.distMM = NaN;  c.nCritical = 0;
        end
        Hh = rheome.differential.helmholtz(J, S);
        e = @(X) sum(X.^2);
        c.solenoidalFrac = e(Hh.Vsol)/max(e(Hh.Vtot), realmin);
        c.irrotationalFrac = e(Hh.Virr)/max(e(Hh.Vtot), realmin);
        c.normalFrac = e(Hh.Vtot - Hh.Virr - Hh.Vsol - Hh.Hresid)/max(e(Hh.Vtot), realmin);
        c.harmFrac = mean(Hh.HarmFrac);           % helmholtz reports this itself
        out.check = c;
    end
end

% Author: Diellor Basha, 2026
