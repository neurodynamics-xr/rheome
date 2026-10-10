classdef tPhaseGradient < matlab.unittest.TestCase
% rheome.flow.phasegradient and rheome.detect.phasesingularity -- the local wavevector of a complex
% scalar field, and the topological charge its phase carries.
%
% The ground truth here is analytic: a planted plane wave has a known k, a planted vortex
% has a known winding number, and a closed surface must carry zero total charge.
%
% Author: Diellor Basha, 2026

    methods (Static)
        function S = flatMesh(n, L)
            % A flat triangulated patch in the z = 0 plane, side L, n x n vertices.
            [X, Y] = meshgrid(linspace(-L/2, L/2, n));
            V = [X(:), Y(:), zeros(numel(X), 1)];
            F = delaunay(V(:,1), V(:,2));
            % delaunay does NOT promise consistent winding; enforce CCW so the patch is a
            % legitimately oriented surface (rheome.detect.phasesingularity requires one).
            d = (V(F(:,2),1)-V(F(:,1),1)).*(V(F(:,3),2)-V(F(:,1),2)) ...
              - (V(F(:,3),1)-V(F(:,1),1)).*(V(F(:,2),2)-V(F(:,1),2));
            F(d < 0, :) = F(d < 0, [1 3 2]);
            S = struct('Vertices', V, 'Faces', F, ...
                       'VertNormals', repmat([0 0 1], size(V,1), 1));
        end
        function S = sphereMesh(n)
            rng(7);
            p = randn(n, 3);  p = p ./ vecnorm(p, 2, 2);
            F = convhull(p(:,1), p(:,2), p(:,3));
            S = struct('Vertices', p, 'Faces', F, 'VertNormals', p);
        end
    end

    methods (Test)

        function planeWaveRecoversItsWavevector(tc)
            % z = exp(i k.x) has grad(phase) = k everywhere. This is the defining case.
            S = tPhaseGradient.flatMesh(40, 0.20);          % 20 cm patch
            k = [2*pi/0.05, 2*pi/0.08, 0];                  % 5 cm and 8 cm wavelengths
            z = exp(1i * (S.Vertices * k(:)));
            out = rheome.flow.phasegradient(z, S);
            kf = squeeze(out.k);                            % [nF x 3]
            tc.verifyEqual(median(kf(:,1)), k(1), 'RelTol', 1e-6);
            tc.verifyEqual(median(kf(:,2)), k(2), 'RelTol', 1e-6);
            tc.verifyLessThan(max(abs(kf(:,3))), 1e-6);     % tangent to the surface
        end

        function wavelengthIsReportedInMillimetres(tc)
            S = tPhaseGradient.flatMesh(40, 0.20);
            lam = 0.05;                                     % 50 mm
            z = exp(1i * (S.Vertices(:,1) * (2*pi/lam)));
            out = rheome.flow.phasegradient(z, S);
            tc.verifyEqual(median(out.wavelength(:)), 50, 'RelTol', 1e-6);
        end

        function manyCyclesAcrossTheMeshNeedNoUnwrapping(tc)
            % THE POINT OF THE LOG-DERIVATIVE FORM. The phase here spans ~40 rad, far beyond
            % one branch of atan2. A route through unwrap() would have to make a global
            % choice; this one is local and never sees a branch cut.
            S = tPhaseGradient.flatMesh(60, 0.20);
            k = 2*pi/0.03;                                  % 30 mm over a 200 mm patch
            z = exp(1i * (S.Vertices(:,1) * k));
            tc.verifyGreaterThan(max(angle(z))-min(angle(z)), 6);       % genuinely wrapped
            out = rheome.flow.phasegradient(z, S);
            kf = squeeze(out.k);
            tc.verifyEqual(median(kf(:,1)), k, 'RelTol', 1e-6);
        end

        function amplitudeGradientSeparatesFromPhaseGradient(tc)
            % grad(log A) and grad(phi) are the real and imaginary halves of one decomposition.
            % Ramp the amplitude along x and the phase along y: neither may leak into the other.
            S = tPhaseGradient.flatMesh(30, 0.10);
            a = 3.0;  ky = 2*pi/0.04;
            z = exp(a * S.Vertices(:,1)) .* exp(1i * ky * S.Vertices(:,2));
            out = rheome.flow.phasegradient(z, S);
            kf = squeeze(out.k);
            tc.verifyLessThan(max(abs(kf(:,1))), 1e-9);          % no amplitude leak into k
            tc.verifyEqual(median(kf(:,2)), ky, 'RelTol', 1e-9);
            ag = squeeze(out.ampgrad);
            tc.verifyEqual(median(ag(:,1)), a, 'RelTol', 1e-9);
            tc.verifyLessThan(max(abs(ag(:,2))), 1e-9);          % no phase leak into amplitude
        end

        function aRealFieldIsRejectedRatherThanSilentlyGivingZero(tc)
            % The phase of a real signal is 0 or pi, so k would come back identically zero and
            % look like a legitimate "no propagation" result. Refuse instead.
            S = tPhaseGradient.flatMesh(10, 0.05);
            tc.verifyError(@() rheome.flow.phasegradient(exp(S.Vertices(:,1)), S), ...
                'flow:phasegradient:real');
        end

        function temporalFrequencyIsRecoveredFromTheSameTrick(tc)
            % omega = Im(dz/dt / z), the time-axis twin of the spatial gradient.
            S = tPhaseGradient.flatMesh(8, 0.05);
            fs = 600;  t = (0:399)/fs;  f0 = 10;
            z = repmat(exp(2i*pi*f0*t), size(S.Vertices,1), 1);
            out = rheome.flow.phasegradient(z, S, 'Rate', fs);
            w = out.omega(:, 3:end-2);                      % interior, away from edge effects
            tc.verifyEqual(median(w(:)), 2*pi*f0, 'RelTol', 1e-3);
        end

        function propagationVelocityCombinesOmegaAndK(tc)
            % A travelling wave exp(i(kx - wt)) moves at w/k along +x.
            S  = tPhaseGradient.flatMesh(30, 0.20);
            fs = 600;  t = (0:299)/fs;  f0 = 10;  lam = 0.06;
            kx = 2*pi/lam;
            z  = exp(1i * (S.Vertices(:,1)*kx - 2*pi*f0*t));
            out = rheome.flow.phasegradient(z, S, 'Rate', fs);
            sp = out.speed(:, 3:end-2);
            tc.verifyEqual(median(sp(:)), lam*f0, 'RelTol', 1e-2);   % 0.6 m/s
            vx = out.velocity(:, 1, 3:end-2);  vy = out.velocity(:, 2, 3:end-2);
            tc.verifyEqual(median(vx(:)), lam*f0, 'RelTol', 1e-2);   % along +x, not -x
            tc.verifyLessThan(abs(median(vy(:))), 1e-6);
        end

        function theEdgeWrappedFormBeatsTheDividedDifferenceAtCoarseResolution(tc)
            % WHY THE ESTIMATOR WAS CHANGED. k = Im(grad z / z) with z taken as the face mean
            % is first order and biased UPWARD, so wavelengths read short and propagation reads
            % fast -- the direction that flatters a travelling-wave claim. Pinned here so a
            % future simplification back to the divided difference fails loudly.
            L = 0.20;  n = 61;
            [X, Y] = meshgrid(linspace(-L/2, L/2, n));
            V = [X(:), Y(:), zeros(numel(X),1)];
            F = delaunay(V(:,1), V(:,2));
            S = struct('Vertices', V, 'Faces', F, 'VertNormals', repmat([0 0 1], size(V,1), 1));
            h = L/(n-1);  fg = rheome.operators.face_gradient(V, F);
            a = F(:,1); b = F(:,2); c = F(:,3);
            for spw = [8 4]
                k = 2*pi/(spw*h);
                z = exp(1i * (V(:,1) * k));
                exact = median(rheome.flow.phasegradient(z, S).kmag);
                zf = (z(a)+z(b)+z(c))/3;
                naive = median(vecnorm([imag((fg.Gx*z)./zf), imag((fg.Gy*z)./zf), ...
                                        imag((fg.Gz*z)./zf)], 2, 2));
                tc.verifyEqual(exact, k, 'RelTol', 1e-9);            % exact at any resolution
                tc.verifyGreaterThan(naive/k - 1, 0);                % and the bias is upward
            end
            % the coarser the mesh the worse it gets: ~3.5% at 8 samples, ~14.6% at 4
            kf = @(spw) 2*pi/(spw*h);
            nb = @(spw) median(vecnorm([imag((fg.Gx*exp(1i*V(:,1)*kf(spw)))./((exp(1i*V(a,1)*kf(spw))+exp(1i*V(b,1)*kf(spw))+exp(1i*V(c,1)*kf(spw)))/3)), ...
                                        imag((fg.Gy*exp(1i*V(:,1)*kf(spw)))./((exp(1i*V(a,1)*kf(spw))+exp(1i*V(b,1)*kf(spw))+exp(1i*V(c,1)*kf(spw)))/3)), ...
                                        imag((fg.Gz*exp(1i*V(:,1)*kf(spw)))./((exp(1i*V(a,1)*kf(spw))+exp(1i*V(b,1)*kf(spw))+exp(1i*V(c,1)*kf(spw)))/3))], 2, 2))/kf(spw) - 1;
            tc.verifyEqual(nb(8), 0.035, 'AbsTol', 0.004);
            tc.verifyEqual(nb(4), 0.146, 'AbsTol', 0.010);
        end

        function aPlantedVortexCarriesUnitCharge(tc)
            % z = (x-x0) + i(y-y0) winds once, counter-clockwise, about (x0,y0).
            S = tPhaseGradient.flatMesh(41, 0.20);
            x0 = 0.013;  y0 = -0.007;                       % off-vertex, off-edge
            z = (S.Vertices(:,1)-x0) + 1i*(S.Vertices(:,2)-y0);
            out = rheome.detect.phasesingularity(z, S);
            tc.verifyEqual(out.total, 1);
            tc.verifyEqual(numel(out.charge), 1);
            tc.verifyEqual(out.charge, 1);
            tc.verifyLessThan(norm(out.pos(1,1:2) - [x0 y0]), 0.012);  % within a triangle
        end

        function conjugatingTheFieldFlipsTheChirality(tc)
            S = tPhaseGradient.flatMesh(41, 0.20);
            z = (S.Vertices(:,1)-0.013) + 1i*(S.Vertices(:,2)+0.007);
            a = rheome.detect.phasesingularity(z, S);
            b = rheome.detect.phasesingularity(conj(z), S);
            tc.verifyEqual(a.total,  1);
            tc.verifyEqual(b.total, -1);
        end

        function aClosedSurfaceCarriesZeroNetCharge(tc)
            % A complex SCALAR field on a closed surface must have total winding 0 -- unlike a
            % vector field, whose indices sum to the Euler characteristic. This is the global
            % consistency check, and it is what distinguishes this from rheome.detect.criticalPoints.
            S = tPhaseGradient.sphereMesh(700);
            rng(11);
            z = complex(randn(size(S.Vertices,1),1), randn(size(S.Vertices,1),1));
            out = rheome.detect.phasesingularity(z, S);
            tc.verifyEqual(out.total, 0);
            tc.verifyGreaterThan(numel(out.charge), 0);     % noise does make singularities
        end

        function netChargeSurvivesVertexNormalsThatLocallyDisagree(tc)
            % REGRESSION. Flipping each face against its averaged vertex normal breaks the
            % edge cancellation that forces sum(charge) = 0, and does it silently. On the real
            % cortex 0.84% of faces disagree (sharp folds) and that was enough to put the net
            % charge off zero in a quarter of frames. Corrupt some normals and demand 0 anyway.
            S = tPhaseGradient.sphereMesh(600);
            rng(3);
            bad = randperm(size(S.Vertices,1), round(0.05*size(S.Vertices,1)));
            S.VertNormals(bad,:) = -S.VertNormals(bad,:);      % 5% point the wrong way
            z = complex(randn(size(S.Vertices,1),1), randn(size(S.Vertices,1),1));
            out = rheome.detect.phasesingularity(z, S);
            tc.verifyEqual(out.total, 0);
        end

        function reversingEveryFaceKeepsTheOutwardChirality(tc)
            % Outward is outward however the faces are wound: reversing all of them flips both
            % the raw winding and the orientation vote, and the two cancel.
            S = tPhaseGradient.flatMesh(41, 0.20);
            z = (S.Vertices(:,1)-0.0123) + 1i*(S.Vertices(:,2)+0.0071);
            a = rheome.detect.phasesingularity(z, S);
            S2 = S;  S2.Faces = S.Faces(:, [1 3 2]);
            b = rheome.detect.phasesingularity(z, S2);
            tc.verifyEqual(a.total, 1);
            tc.verifyEqual(b.total, a.total);
        end

        function anExactlyPiPhaseStepStillCancelsAcrossItsTwoFaces(tc)
            % REGRESSION on the wrap. A core placed exactly antipodal to two vertices makes one
            % edge step exactly pi; with mod()-based wrapping that edge stopped cancelling and
            % the planted vortex read as charge 0. Grid 5 mm, core on a half-step diagonal.
            S = tPhaseGradient.flatMesh(41, 0.20);
            z = (S.Vertices(:,1)-0.013) + 1i*(S.Vertices(:,2)-0.007);
            phi = angle(z);  F = S.Faces;
            wr = @(x) atan2(sin(x), cos(x));
            steps = abs(wr(phi(F(:,2))-phi(F(:,1))));
            tc.verifyGreaterThan(max(steps), pi - 1e-9);        % the degenerate step exists
            out = rheome.detect.phasesingularity(z, S);
            tc.verifyEqual(out.total, 1);                       % and it is still found
        end

        function asmoothFieldWithNoVortexHasNoSingularities(tc)
            S = tPhaseGradient.flatMesh(40, 0.20);
            z = exp(1i * (S.Vertices(:,1) * (2*pi/0.08)));
            out = rheome.detect.phasesingularity(z, S);
            tc.verifyEmpty(out.charge);
            tc.verifyEqual(out.total, 0);
        end

        function chargeIsInvariantToAmplitudeScaling(tc)
            % The reason to prefer winding over raw curl: multiplying the field by any
            % positive amplitude profile cannot move or destroy a topological charge.
            S = tPhaseGradient.flatMesh(41, 0.20);
            z = (S.Vertices(:,1)-0.013) + 1i*(S.Vertices(:,2)-0.007);
            A = 0.01 + exp(-sum(S.Vertices(:,1:2).^2,2)/0.002);   % strong, positive, varying
            a = rheome.detect.phasesingularity(z, S);
            b = rheome.detect.phasesingularity(A .* z, S);
            tc.verifyEqual(b.total, a.total);
            tc.verifyEqual(b.pos, a.pos, 'AbsTol', 1e-12);
        end

    end
end
