classdef tCorticalScales < matlab.unittest.TestCase
% The spatial pyramid on a surface: a nested partition (rheome.geom.tree) and a scale ruler
% (rheome.detect.blobscale).
%
% These are the two halves the time axis already has, and they divide the same way. The tree
% is disjoint and nested, so sums over its vertices roll up exactly -- it is bookkeeping. The
% scale selection is smooth and overlapping, so nothing merges -- it is measurement. What is
% asserted is that the partition really partitions, that its levels scale the way a
% two-dimensional bisection must, and that the ruler returns the size that was planted.
%
% Author: Diellor Basha, 2026

    properties
        S; L; M; basis; R = 0.07            % a head-sized sphere, so wavelengths are in mm
    end

    methods (TestClassSetup)
        function mesh(tc)
            [V, F] = rheome.geom.icosphere(3);                         % 642 vertices
            tc.S = struct('Vertices', V * tc.R, 'Faces', F, 'nV', size(V,1), 'nF', size(F,1));
            [tc.L, tc.M] = rheome.operators.laplace_beltrami(tc.S.Vertices, tc.S.Faces);
            tc.basis = rheome.eigen.modes(tc.L, tc.M, 120);
        end
    end

    methods (Test)

        % ---------------- the partition ----------------

        function theLeavesPartitionTheSurfaceExactly(tc)
            % ⭐ The property the whole bookkeeping side rests on: disjoint and covering, so a
            % parent's sum is its children's sums.
            T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=4);
            lv = T(T.is_leaf, :);
            all_ = sort([lv.members{:}]);
            tc.verifyEqual(all_, 1:tc.S.nV);
            tc.verifyEqual(numel(unique(all_)), tc.S.nV);         % disjoint, not merely covering
            for i = find(~T.is_leaf)'
                kids = T(T.parent_id == T.node_id(i), :);
                tc.verifyEqual(height(kids), 2);
                tc.verifyEqual(sort([kids.members{:}]), sort(T.members{i}));
            end
        end

        function areaAddsUpTheTree(tc)
            T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=4);
            total = full(sum(sum(tc.M)));
            tc.verifyEqual(T.area(T.depth == 0), total, 'RelTol', 1e-10);
            tc.verifyEqual(sum(T.area(T.is_leaf)), total, 'RelTol', 1e-10);
            for d = 1:max(T.depth)
                tc.verifyEqual(sum(T.area(T.depth == d)), total, 'RelTol', 1e-10, sprintf('depth %d', d));
            end
        end

        function eachLevelHalvesTheAreaAndTheDiameterFallsBySqrtTwo(tc)
            % ⚠ A BISECTION IS NOT AN OCTAVE. Halving the AREA divides the diameter by
            % sqrt(2), so a dyadic ladder in wavelength is every OTHER level -- which is the
            % factor of four per octave Weyl's law predicts for a two-dimensional domain.
            T = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=4);
            for d = 1:max(T.depth)
                a0 = mean(T.area(T.depth == d-1));  a1 = mean(T.area(T.depth == d));
                tc.verifyEqual(a1 / a0, 0.5, 'RelTol', 0.05, sprintf('area at depth %d', d));
                d0 = mean(T.diameter(T.depth == d-1));  d1 = mean(T.diameter(T.depth == d));
                tc.verifyEqual(d1 / d0, 1/sqrt(2), 'RelTol', 0.05, sprintf('diameter at depth %d', d));
            end
        end

        function disconnectedPartsAreSeparatedBeforeAnyEigensolve(tc)
            % ⚠ Two hemispheres make the Fiedler value zero and the eigenproblem degenerate.
            % A node spanning components is split by component instead.
            [V, F] = rheome.geom.icosphere(2);
            V2 = [V * tc.R; V * tc.R + [0.5 0 0]];
            F2 = [F; F + size(V,1)];
            S2 = struct('Vertices', V2, 'Faces', F2, 'nV', size(V2,1), 'nF', size(F2,1));
            T = rheome.geom.tree(S2, MaxDepth=2);
            kids = T(T.parent_id == 1, :);
            tc.verifyEqual(height(kids), 2);
            for i = 1:2
                m = kids.members{i};
                tc.verifyTrue(all(m <= size(V,1)) || all(m > size(V,1)), 'a child straddles the two spheres');
            end
        end

        function theStoppingRulesStop(tc)
            T2 = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MaxDepth=2);
            tc.verifyEqual(max(T2.depth), 2);
            tc.verifyEqual(height(T2), 7);                        % 1 + 2 + 4
            total = full(sum(sum(tc.M)));
            T3 = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MinArea=total/5);
            tc.verifyTrue(all(T3.area(~T3.is_leaf) > total/5));
            T4 = rheome.geom.tree(tc.S, L=tc.L, M=tc.M, MinVertices=200);
            tc.verifyTrue(all(T4.n_vertices(~T4.is_leaf) > 200));
        end

        % ---------------- the ruler ----------------

        function thePlantedSizeIsTheSizeItReports(tc)
            % ⭐ The scale-selection claim, against a blob of known aperture.
            for wl = [0.05 0.09]
                t = (wl/(2*pi))^2;
                d = zeros(tc.S.nV, 1);  d(1) = 1;
                x = tc.basis.Phi * (exp(-tc.basis.Lambda * t) .* (tc.basis.Phi' * (tc.M * d)));
                out = rheome.detect.blobscale(x, tc.basis, Wavelengths=[0.02 0.3], NumScales=16);
                step = out.wavelengths(2) / out.wavelengths(1);
                tc.verifyEqual(out.wavelength(1) / wl, 1, 'RelTol', step - 1 + 0.05, ...
                    sprintf('planted %.0f mm, read %.0f mm', 1e3*wl, 1e3*out.wavelength(1)));
            end
        end

        function theLadderAndItsDerivedQuantitiesAgree(tc)
            x = randn(tc.S.nV, 3);
            out = rheome.detect.blobscale(x, tc.basis, Wavelengths=[0.02 0.2], NumScales=9);
            tc.verifyNumElements(out.scales, 9);
            tc.verifyEqual(out.wavelengths, 2*pi*sqrt(out.scales), 'RelTol', 1e-12);
            tc.verifyEqual(out.aperture, sqrt(2*out.scale), 'RelTol', 1e-12);
            tc.verifySize(out.scale, [tc.S.nV, 3]);
            tc.verifyTrue(all(ismember(out.wavelength(:), out.wavelengths)));
            r = out.wavelengths(2:end) ./ out.wavelengths(1:end-1);
            tc.verifyEqual(r, repmat(r(1), size(r)), 'RelTol', 1e-9);   % log-spaced
        end

        function theResolutionFloorIsHonoured(tc)
            % ⚠ Below the instrument's aperture the selected scale is the estimator's point
            % spread, not the cortex, so a floor removes those rungs rather than reporting them.
            x = randn(tc.S.nV, 1);
            out = rheome.detect.blobscale(x, tc.basis, Wavelengths=[0.02 0.3], NumScales=12, MinWavelength=0.1);
            tc.verifyGreaterThanOrEqual(min(out.wavelengths), 0.1);
            tc.verifyGreaterThanOrEqual(min(out.wavelength(:)), 0.1);
            tc.verifyError(@() rheome.detect.blobscale(x, tc.basis, Wavelengths=[0.02 0.05], MinWavelength=0.5), ...
                           'detect:blobscale:floor');
        end

        function aNegativeEigenvalueIsClampedAndSaidSo(tc)
            % ⚠ THE TRAP: exp(-lambda*t) with lambda < 0 GROWS, so one spurious mode pins
            % every selection to the coarsest rung. Measured on the real cortex: eigs with
            % 'smallestabs' returned a -1.6e4 mode and every blob read as 400 mm.
            b = tc.basis;  b.Lambda(5) = -1e4;
            x = randn(tc.S.nV, 1);
            tc.verifyWarning(@() rheome.detect.blobscale(x, b, Wavelengths=[0.02 0.2]), 'detect:blobscale:negative');
            w = warning('off', 'detect:blobscale:negative');
            out = rheome.detect.blobscale(x, b, Wavelengths=[0.02 0.2]);
            warning(w);
            tc.verifyTrue(all(isfinite(out.response(:))));
        end

        function itWorksOnAComplexFieldAndWithAnotherKernel(tc)
            z = complex(randn(tc.S.nV, 1), randn(tc.S.nV, 1));
            out = rheome.detect.blobscale(z, tc.basis, Wavelengths=[0.02 0.2], NumScales=8);
            tc.verifyTrue(isreal(out.response));
            tc.verifyTrue(all(out.response >= 0));
            og = rheome.detect.blobscale(z, tc.basis, Wavelengths=[0.02 0.2], NumScales=8, Kernel="diffgauss");
            tc.verifySize(og.scale, size(out.scale));
            oh = rheome.detect.blobscale(z, tc.basis, Wavelengths=[0.02 0.2], NumScales=8, ...
                                  Kernel=@(l, t) (t*l) .* exp(-t*l), Gamma=0);
            tc.verifyEqual(oh.wavelength, out.wavelength);           % the same kernel by hand
            tc.verifyError(@() rheome.detect.blobscale(z, tc.basis, Kernel="ricker"), 'detect:blobscale:kernel');
            tc.verifyError(@() rheome.detect.blobscale(randn(3,1), tc.basis), 'detect:blobscale:size');
        end

    end
end

% Author: Diellor Basha, 2026
