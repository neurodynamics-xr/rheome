classdef tChargeDensity < matlab.unittest.TestCase
% rheome.detect.chargedensity -- topological charge read at a chosen spatial aperture.
%
% Ground truth: on a sphere the field x+iy has exactly two zeros, one at each pole, carrying
% +1 and -1. A correct aperture-integrated charge must find them, hold their value across a
% range of apertures, and integrate to zero over the closed surface.
%
% Author: Diellor Basha, 2026

    methods (Static)
        function [S, basis] = sphereBasis(nsub, R, K)
            % A REGULAR icosphere, not random points. A convhull of random points gives a mesh
            % whose local density varies several-fold, and the eigenbasis inherits that: the
            % same planted core read very differently at the two poles (0.63 vs -0.77, with a
            % sign flip at the largest aperture) purely from mesh irregularity. Subdivision
            % keeps every triangle within a few percent of the same size.
            t = (1+sqrt(5))/2;
            V = [-1 t 0; 1 t 0; -1 -t 0; 1 -t 0; 0 -1 t; 0 1 t; 0 -1 -t; 0 1 -t;
                  t 0 -1; t 0 1; -t 0 -1; -t 0 1];
            F = [1 12 6; 1 6 2; 1 2 8; 1 8 11; 1 11 12; 2 6 10; 6 12 5; 12 11 3; 11 8 7;
                 8 2 9; 4 10 5; 4 5 3; 4 3 7; 4 7 9; 4 9 10; 5 10 6; 3 5 12; 7 3 11;
                 9 7 8; 10 9 2];
            V = V ./ vecnorm(V,2,2);
            for k = 1:nsub
                nV = size(V,1);  E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
                [Eu,~,ic] = unique(sort(E,2), 'rows');
                mid = nV + ic;
                V = [V; (V(Eu(:,1),:) + V(Eu(:,2),:))/2];
                m1 = mid(1:size(F,1)); m2 = mid(size(F,1)+1:2*size(F,1));
                m3 = mid(2*size(F,1)+1:end);
                F = [F(:,1) m1 m3; F(:,2) m2 m1; F(:,3) m3 m2; m1 m2 m3];
                V = V ./ vecnorm(V,2,2);
            end
            V = R * V;
            S = struct('Vertices', V, 'Faces', F, 'VertNormals', V ./ vecnorm(V,2,2));
            % LUMPED mass, so M is diagonal and the column normalisation below really is
            % mass-orthonormalisation. With the default Galerkin M the normalisation ignores
            % the off-diagonal terms, Phi'*M*Phi ~= I, and Phi*E*Phi' STOPS BEING THE HEAT
            % KERNEL -- two cores 5.3 mm from their nearest vertex then read 0.29 and 5.10.
            [L, M] = rheome.operators.laplace_beltrami(V, F, 'lumped');
            [Phi, D] = eigs(L, M, K, 'smallestabs');
            w = full(diag(M));
            Phi = Phi ./ sqrt(sum(w .* Phi.^2, 1));
            basis = struct('Phi', Phi, 'Lambda', diag(D));
        end
        function v = nearestTo(S, out, want)
            % The vertex closest to the detected core of the requested sign. An icosphere has
            % NO vertex at the z-pole -- the nearest is 32 deg away, ~33 mm on a 60 mm sphere --
            % so reading q at max(z) samples a core from more than an aperture away.
            f = find(out.perFace(:,1) == want);
            F = S.Faces;  V = S.Vertices;
            cen = (V(F(f,1),:) + V(F(f,2),:) + V(F(f,3),:)) / 3;
            [~, v] = min(vecnorm(V - cen(1,:), 2, 2));
        end
    end

    methods (Test)

        function xPlusIYOnASphereHasExactlyTwoOppositeCores(tc)
            [S, basis] = tChargeDensity.sphereBasis(3, 0.06, 250);
            z = S.Vertices(:,1) + 1i*S.Vertices(:,2);
            out = rheome.detect.chargedensity(z, S, basis, 26);
            tc.verifyEqual(sort(out.perFace(out.perFace(:,1)~=0, 1)).', [-1 1]);
            tc.verifyEqual(sum(out.perFace(:,1)), 0);
        end

        function aUnitCoreReadsAsUnitChargeNotAsADensity(tc)
            % REGRESSION. Phi*E*Phi'*c is a DENSITY: at a unit core it reads ~1/(4*pi*tau) --
            % thousands -- and GROWS as the aperture shrinks, so a fixed +-0.5 threshold passed
            % 797 of 900 vertices on pure noise. Normalising by the kernel's own centre value
            % K(v,v) makes an enclosed unit charge read order 1 whatever the aperture.
            [S, basis] = tChargeDensity.sphereBasis(3, 0.06, 250);
            z = S.Vertices(:,1) + 1i*S.Vertices(:,2);
            out = rheome.detect.chargedensity(z, S, basis, [26 32]);
            vp = tChargeDensity.nearestTo(S, out, +1);
            vm = tChargeDensity.nearestTo(S, out, -1);
            tc.verifyGreaterThan(out.qmed(vp,1),  0.4);
            tc.verifyLessThan(   out.qmed(vm,1), -0.4);
            tc.verifyLessThan(max(abs(out.q(:))), 20);
        end

        function aGenuineCoreIsFlatterAcrossAperturesThanNoise(tc)
            % THE DISCRIMINATOR. An isolated unit core encloses ~1 at every aperture below its
            % isolation scale and HOLDS; the dense mixed-sign cores of a random field never
            % reach 1 and do not hold. Measured at r = 26/32/40 mm on a 60 mm sphere: core
            % 0.894/0.929/0.954 flatness 0.937, noise max 0.96/0.93/0.78 flatness 0.528.
            [S, basis] = tChargeDensity.sphereBasis(3, 0.06, 250);
            radii = [26 32 40];
            zc = S.Vertices(:,1) + 1i*S.Vertices(:,2);
            oc = rheome.detect.chargedensity(zc, S, basis, radii);
            vp = tChargeDensity.nearestTo(S, oc, +1);

            rng(21);
            zr = complex(randn(size(S.Vertices,1),1), randn(size(S.Vertices,1),1));
            orr = rheome.detect.chargedensity(zr, S, basis, radii);
            tc.verifyGreaterThan(orr.nCore(1), 10);              % noise DOES make raw cores
            tc.verifyGreaterThan(oc.qmed(vp,1), 0.8);            % the real one reads ~1
            % flatness is bounded by 1, so compare by margin, not by ratio
            tc.verifyGreaterThan(oc.flatness(vp,1), 0.9);
            tc.verifyLessThan(median(orr.flatness(:,1)), 0.65);
            % Neither half separates alone at the tail -- the background's best vertex matches
            % the core on EITHER axis. Only the joint criterion is selective, so pin the rate
            % over several realisations rather than one (a single draw is noisy at ~0.6%).
            tc.verifyTrue(abs(oc.qmed(vp,1)) >= 0.8 && oc.flatness(vp,1) >= 0.90);
            hit = 0;  tot = 0;
            for k = 1:5
                rng(100+k);
                zk = complex(randn(size(S.Vertices,1),1), randn(size(S.Vertices,1),1));
                ok = rheome.detect.chargedensity(zk, S, basis, radii);
                hit = hit + sum(abs(ok.qmed(:,1)) >= 0.8 & ok.flatness(:,1) >= 0.90);
                tot = tot + size(S.Vertices,1);
            end
            tc.verifyLessThan(hit/tot, 0.02);                    % measured ~0.6% per vertex
        end

        function anInconsistentlyOrientedMeshIsRefused(tc)
            [S, basis] = tChargeDensity.sphereBasis(2, 0.06, 100);
            S.Faces(1,:) = S.Faces(1,[1 3 2]);
            z = complex(randn(size(S.Vertices,1),1), randn(size(S.Vertices,1),1));
            tc.verifyError(@() rheome.detect.chargedensity(z, S, basis, 26), ...
                'detect:chargedensity:orientation');
        end

        function theApertureFloorIsReportedNotSilentlyImposed(tc)
            % Below the finest wavelength the basis represents, the requested radius is NOT what
            % you get -- at 250 modes on this sphere the floor is ~20 mm. .resolved must fall
            % away so the caller can see it rather than trust a number the basis cannot express.
            [S, basis] = tChargeDensity.sphereBasis(3, 0.06, 250);
            z = S.Vertices(:,1) + 1i*S.Vertices(:,2);
            out = rheome.detect.chargedensity(z, S, basis, [3 12 45]);
            tc.verifyLessThan(out.resolved(1), out.resolved(3));
            tc.verifyGreaterThan(out.resolved(3), 0.9);
            tc.verifyLessThan(out.resolved(1), 0.5);
        end

    end
end
