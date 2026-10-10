function out = directionfield(name, varargin)
% FLOW.DIRECTIONFIELD  Plant singularities, solve a trivial connection, integrate, rotate, lift.
%
%   out = rheome.flow.directionfield('sub01')                        % two +1 poles
%   out = rheome.flow.directionfield(name, Faces=[f1 f2 f3 f4], Charges=[1 -1 1 1])   % a +1/-1 pair, plus chi
%   out = rheome.flow.directionfield(name, Rotate=0)                         % the un-rotated field
%
% ⭐⭐ THE DESIGN ROUTE, which is what geometry-processing-js / Geometry Central do. You do NOT derive
% a field and then detect what came out; you PRESCRIBE the singularities and solve for the field
% that has them:
%
%   1. choose faces and integer indices k_f, summing to chi        (the singularities)
%   2. solve the TRIVIAL CONNECTION (Crane, Desbrun & Schroeder 2010), rheome.operators.trivial_connection:
%          least-norm 1-form x with  d1*x = 2*pi*k - K,  K_f = wrap(d1*argR)_f  the face curvature
%   3. INTEGRATE the adjusted transport over a spanning tree  ->  phi(v), a unit direction field
%   4. ROTATE by Rotate (default 90 degrees): z = exp(i*(phi + Rotate))
%   5. LIFT to 3 components through the connection's own frames:  J = Re(z)*e1 + Im(z)*e2
%
% ⚠⚠ THE INDICES SUM TO chi, AND AN EARLIER VERSION SAID ZERO. It solved d1*x = 2*pi*k - d1*argR
% without the wrap; d1*argR = K + 2*pi*n carries integer chart parts n_f (sum -chi), so the budget
% looked like zero and the field's index was k - n: a "+1/-1 pair" wound on 2560 faces. Its checks
% passed because they measured the ADJUSTED connection's holonomy, chart parts included, which is
% exactly what was prescribed; the field's own winding (.check.winding*) said otherwise, and a test
% pinned that as a detector limitation. It was not: with the wrap the per-face winding of the field
% is nonzero at exactly the prescribed faces, so detection and prescription now agree.
%
% ⭐ SO A VORTEX/ANTIVORTEX PAIR CAN BE PLANTED -- next to the chi budget, which has to go somewhere:
% Faces=[a b p q], Charges=[1 -1 1 1] puts the pair at a, b and two +1 poles at p, q. The default,
% two +1 faces far apart, is the pole frame of rheome.operators.gauge(Method="trivial").
%
% ⚠ THE UN-ROTATED FIELD IS NOT LITERALLY IRROTATIONAL, and calling it that is loose. A parallel unit
% field has zero covariant derivative, so it is as close to constant as the surface allows -- its
% singularities are the prescribed ones and its curl is whatever the curvature forces. Rotating by 90
% degrees exchanges the roles of div and curl locally; it does not turn an exactly-curl-free field into
% an exactly-divergence-free one, because neither is exact on a curved discrete surface.
%
% ⚠ A UNIT FIELD DOES NOT DECAY. |z| = 1 everywhere, so this covers the whole hemisphere. Pass Envelope
% to localise it; an envelope that is radial about a singularity adds no divergence, because grad(A) is
% then perpendicular to the rotated field.
%
% NAME-VALUE
%   Faces     [1 x n] face indices for the singularities; default two well-separated faces
%   Charges   [1 x n] integers summing to chi (2 on a closed hemisphere); default [1 1]
%   Rotate    radians added to phi; pi/2 (default) = the rotated/vortex field, 0 = the parallel one
%   Envelope  "none" (default) | "radial"   -- radial decays about the FIRST singularity
%   EnvelopeMM  the radial envelope's scale, default 140
%   Hemi ("L")   Bases   Check (true)
%
% OUTPUT (struct out)
%   .J [3nV x 1]  .z .phi  .faces .charges  .check: holonomyAtPrescribed, holonomyElsewhereMax
%   (the adjusted connection's index per face, in turns, chart parts removed), transportResidual,
%   windingAtPrescribed, windingSumAll, windingNonzeroCount (the FIELD's own per-face winding),
%   tangentialFrac, curlOverDiv
%
% See also: rheome.operators.trivial_connection, rheome.operators.connection_laplacian, rheome.operators.gauge, rheome.flow.seedvortex, rheome.flow.seedrotor,
%           rheome.detect.criticalPoints
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Faces', [], @isnumeric);
    p.addParameter('Charges', [], @isnumeric);
    p.addParameter('Rotate', pi/2, @isscalar);
    p.addParameter('Envelope', "none");
    p.addParameter('EnvelopeMM', 140, @isscalar);
    p.addParameter('Hemi', "L");
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Connection', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Check', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    if isempty(o.Bases), B = rheome.load.bases(name); else, B = o.Bases; end
    H = B.(char(o.Hemi));  S = H.S;  V = S.Vertices;  F = double(S.Faces);
    nV = size(V,1);  nF = size(F,1);
    if isempty(o.Connection), C = rheome.operators.connection_laplacian(V, F); else, C = o.Connection; end
    wr = @(x) mod(x+pi, 2*pi) - pi;

    %% the prescribed singularities
    fs = o.Faces;  kk = o.Charges;
    if isempty(fs)
        fc = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:))/3;
        [~, i1] = max(vecnorm(fc - mean(fc,1), 2, 2));
        [~, i2] = max(vecnorm(fc - fc(i1,:), 2, 2));
        fs = [i1 i2];
    end
    nE = size(unique(sort([F(:,[2 3]); F(:,[3 1]); F(:,[1 2])], 2), 'rows'), 1);
    chi = nV - nE + nF;
    if isempty(kk), kk = ones(1, numel(fs)); kk(end) = chi - (numel(fs)-1); end
    if sum(kk) ~= chi
        error('flow:directionfield:charges', ...
            ['Charges must sum to chi = %d (Poincare-Hopf), not %d: the face curvatures sum to ' ...
             '2*pi*chi, so another total is inconsistent.'], chi, sum(kk));
    end

    %% solve the trivial connection, then integrate it
    T = rheome.operators.trivial_connection(V, F, C, fs, kk);
    phi = T.phase;  Ax = T.Ax;  d1 = T.d1;  argE2 = T.argE + T.x;

    %% rotate, envelope, lift
    z = exp(1i*(phi + o.Rotate));
    if lower(string(o.Envelope)) == "radial"
        [~, vs] = min(vecnorm(V - mean(V(F(fs(1),:),:),1), 2, 2));
        d = rheome.geom.geodesic(S, vs);
        s = o.EnvelopeMM*1e-3/(2*pi);
        z = z .* ((d/s) .* exp(-(d.^2)/(2*s^2)));
    end
    JJ = real(z).*C.e1 + imag(z).*C.e2;
    JJ = JJ / max(vecnorm(JJ, 2, 2));
    J = reshape(JJ', [], 1);

    out = struct('J', J, 'z', z, 'phi', phi, 'faces', fs, 'charges', kk, ...
                 'rotate', o.Rotate, 'check', struct());

    %% verify against the PRESCRIPTION, not with a detector
    if o.Check
        c = struct();
        hol = (d1*argE2)/(2*pi) - T.n;                     % index per face, in turns, chart parts removed
        c.holonomyAtPrescribed  = hol(fs)';
        other = true(nF,1);  other(fs) = false;
        c.holonomyElsewhereMax  = max(abs(hol(other)));
        [ii, jj] = find(triu(C.A ~= 0, 1));
        res = abs(wr(phi(jj) - phi(ii) - full(Ax(sub2ind([nV nV], ii, jj)))));
        c.transportResidual = max(res);
        % the field's own winding at the prescribed faces, from principal values
        a = F(:,1); b = F(:,2); cc = F(:,3);
        aR = @(r1,c1) angle(full(C.Rt(sub2ind([nV nV], r1, c1))));
        q = ( wr(phi(b)-phi(a)-aR(a,b)) + wr(phi(cc)-phi(b)-aR(b,cc)) ...
            + wr(phi(a)-phi(cc)-aR(cc,a)) ) / (2*pi);
        c.windingAtPrescribed = round(q(fs))';
        c.windingSumAll = sum(round(q));
        c.windingNonzeroCount = nnz(round(q));
        c.tangentialFrac = 1 - sum(sum(JJ.*C.normal,2).^2)/max(sum(JJ(:).^2), realmin);
        fg = rheome.operators.face_gradient(V, F);
        av = full(sum(H.lbo.Mass, 2));  nrm = @(y) sqrt(sum(av.*y.^2));
        c.curlOverDiv = nrm(rheome.differential.curl(J,S,fg)) / max(nrm(rheome.differential.divergence(J,S,fg)), realmin);
        out.check = c;
    end
end

% Author: Diellor Basha, 2026
