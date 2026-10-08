function out = phasegradient(z, S, varargin)
% FLOW.PHASEGRADIENT  Local wavevector, frequency and propagation velocity of a complex field.
%
%   out = rheome.flow.phasegradient(z, S)
%   out = rheome.flow.phasegradient(z, S, 'Rate', fs)
%
% Takes the ANALYTIC field z(v,t) -- the complex curl, or any complex scalar on the mesh --
% and returns where its phase is heading, in space and in time. Writing z = A*exp(i*phi):
% k = grad(phi) is the local wavevector, |k| = 2*pi/wavelength, and with omega = d(phi)/dt
% the pattern's timing sweeps past a vertex at velocity omega*k/|k|^2.
%
% ⭐ EDGE-WRAPPED DIFFERENCES, NOT A DIVIDED DIFFERENCE OF z. Along an edge e the wrapped
% phase step wrap(phi_j - phi_i) equals k.e EXACTLY for a plane wave, so solving the two
% independent edge steps of a triangle for the in-plane k is exact -- to machine precision,
% at any resolution, for any number of cycles across the mesh (tPhaseGradient plants a wave
% spanning 40 rad and recovers k to 1e-12).
%
% The obvious alternative -- k = Im(grad z / z) with grad z from rheome.operators.face_gradient --
% is only FIRST ORDER, and biased in a direction that flatters the result: with z on a face
% taken as the mean of its vertices it returns (2/h)*tan(k*h/2), an OVERestimate of |k|, so
% wavelengths read SHORT and apparent propagation reads FAST. Measured on a flat mesh, in
% samples per wavelength: 20 -> +0.6%, 10 -> +2.2%, 8 -> +3.5%, 6 -> +6.3%, 4 -> +14.6%,
% 3 -> +24.1%, against 0.00% for the edge-wrapped form at every one. tPhaseGradient pins it.
%
% ⚠ THE ALIASING LIMIT IS REAL AND IT IS PER EDGE. Wrapping cannot distinguish k.e from
% k.e + 2*pi, so an edge stepping more than pi is indistinguishable from one stepping
% backwards. .aliased reports the fraction of faces with any edge step past 0.9*pi -- read it
% before believing a short wavelength. On a 3 mm cortical mesh this binds below ~6 mm.
%
% k IS PER FACE (a triangle fixes one in-plane vector); omega is per VERTEX (time is a
% per-vertex axis). .speed and .velocity live on faces, with omega averaged onto them.
%
% WHERE THE DESCRIPTION ITSELF BREAKS DOWN. k is a local wavevector only where phi varies
% smoothly. At a phase singularity it does not, and |k| diverges as the core is approached.
% That divergence is the detector, not a defect -- see rheome.detect.phasesingularity.
%
% INPUTS:
%   z  [nV x nT] complex analytic field on the mesh vertices
%   S  surface (.Vertices, .Faces)
%   'Rate'          Hz, required for .omega/.speed/.velocity
%   'MinAmplitude'  relative to the median |z|; below it phase is meaningless -> NaN.
%                   Default 0 (off). Set it for real data, where |z| does reach 0.
%
% OUTPUT (struct out):
%   .k          [nF x 3 x nT]  local wavevector, rad/m, tangent to the face
%   .kmag       [nF x nT]      |k|
%   .wavelength [nF x nT]      2*pi/|k| in MM -- travel this far to be a cycle out of step
%   .ampgrad    [nF x 3 x nT]  grad(log A), 1/m -- the real part of the same decomposition
%   .omega      [nV x nT]      rad/s        (Rate given)
%   .speed      [nF x nT]      |omega|/|k|, m/s
%   .velocity   [nF x 3 x nT]  omega*k/|k|^2
%   .aliased    [1 x nT]       fraction of faces with an edge step past 0.9*pi
%
% See also: rheome.detect.phasesingularity, rheome.operators.face_gradient, rheome.flow.curl
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Rate', []);
    p.addParameter('MinAmplitude', 0);
    p.parse(varargin{:});
    o = p.Results;

    if isreal(z)
        error('flow:phasegradient:real', ...
            ['z is real. This needs the ANALYTIC field -- apply hilbert() in time, or take ' ...
             'the complex coefficients from rheome.flowpage (both Methods already return them).']);
    end

    V = double(S.Vertices);  F = double(S.Faces);
    a = F(:,1);  b = F(:,2);  c = F(:,3);
    nF = size(F,1);  nT = size(z,2);

    e1 = V(b,:) - V(a,:);                                  % [nF x 3] two spanning edges
    e2 = V(c,:) - V(a,:);
    phi = angle(z);
    wr  = @(x) atan2(sin(x), cos(x));      % antisymmetric -- see below
    d1  = wr(phi(b,:) - phi(a,:));                         % [nF x nT] exact k.e1 for a plane wave
    d2  = wr(phi(c,:) - phi(a,:));

    % Solve [e1; e2] k = [d1; d2] per face for the in-plane k. Closed form via the 2x2 Gram,
    % which is the normal equation of the 2x3 system and stays in the face plane by
    % construction -- no explicit tangent frame, no normal component to discard.
    g11 = sum(e1.*e1, 2);  g12 = sum(e1.*e2, 2);  g22 = sum(e2.*e2, 2);
    det = max(g11.*g22 - g12.^2, realmin);
    c1  = ( g22.*d1 - g12.*d2) ./ det;                     % coefficients on e1, e2
    c2  = (-g12.*d1 + g11.*d2) ./ det;
    k   = permute(cat(3, e1, e2), [1 3 2]);                % [nF x 2 x 3]
    kk  = zeros(nF, 3, nT);
    for d = 1:3
        kk(:,d,:) = c1 .* e1(:,d) + c2 .* e2(:,d);
    end

    % Amplitude: grad(log A) is exact for an exponential profile and needs no wrapping.
    fg = rheome.operators.face_gradient(V, F);
    la = log(max(abs(z), realmin));
    ag = permute(cat(3, fg.Gx*la, fg.Gy*la, fg.Gz*la), [1 3 2]);

    zf = (z(a,:) + z(b,:) + z(c,:)) / 3;
    if o.MinAmplitude > 0
        bad = abs(zf) < o.MinAmplitude * median(abs(zf(:)), 'omitnan');
        kk(repmat(permute(bad, [1 3 2]), 1, 3, 1)) = NaN;
        ag(repmat(permute(bad, [1 3 2]), 1, 3, 1)) = NaN;
    end

    out.k          = kk;
    out.ampgrad    = ag;
    out.kmag       = reshape(vecnorm(kk, 2, 2), nF, nT);
    out.wavelength = 2*pi ./ out.kmag * 1000;
    d3 = wr(phi(c,:) - phi(b,:));
    out.aliased = mean(max(abs(cat(3,d1,d2,d3)), [], 3) > 0.9*pi, 1);

    if ~isempty(o.Rate)
        pw = wr(diff(phi, 1, 2));                          % wrapped step between samples
        w  = [pw(:,1), (pw(:,1:end-1) + pw(:,2:end))/2, pw(:,end)] * o.Rate;
        out.omega    = w;
        wf = (w(a,:) + w(b,:) + w(c,:)) / 3;
        out.speed    = abs(wf) ./ out.kmag;
        out.velocity = kk .* permute(wf ./ out.kmag.^2, [1 3 2]);
    end
    out.nFace = nF;  out.nVert = size(z,1);
end

% Author: Diellor Basha, 2026
