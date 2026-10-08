classdef tScatter < matlab.unittest.TestCase
% rheome.filters.scatter -- graph wavelet scattering, |W| cascaded through the frame and pooled.
%
% The properties asserted here are the ones that are PROVABLE, not merely plausible:
% 1-homogeneity, non-negativity, the coarser-second ordering rule, and block invariance
% (the memory bound must not change the answer). Plus two physical checks: a constant
% field must fall through every band-pass member, and a planted atom of known width must
% peak on the member matched to it.
%
% Author: Diellor Basha, 2026

    properties (Constant)
        MB = 5e8;     % the production budget, so the shape tests exercise one block
    end

    methods (Static)
        function fx = fixture()
            % ico4 (2562 vertices), K = 150 on a 100 mm sphere. Deliberately small: peak
            % memory is ~86 kB per time column, so the suite can never pressure the machine.
            % K = 150 is the largest that does NOT trip rheome.eigen.modes' verts-per-half-wave
            % warning here (3.5 v/hw); K = 200 does. ico3/K = 60 was too coarse to test a
            % scale readout at all -- its usable band was two members wide.
            persistent c
            if ~isempty(c), fx = c; return; end
            R = 0.100;
            [V, F] = rheome.geom.icosphere(4);
            V = R * (V ./ vecnorm(V, 2, 2));
            [L, Mm] = rheome.operators.laplace_beltrami(V, F, 'galerkin');
            b = rheome.eigen.modes(L, Mm, 150);
            fr = rheome.filters.frame('mexhat', [], b.Lambda);   % [] -> derived density; Nf=4 warns
            fx = struct('basis', b, 'V', V, 'F', F, 'frame', fr, ...
                        'M', numel(fr.g), 'nV', size(V,1), 'R', R);
            c = fx;
        end

        function t = tol(ref, k)
            % ⚠ NEVER RelTol ALONE HERE. S2 entries run from O(1) down to 1e-19, and BLAS
            % sums a 1-column block in a different order from a 17-column one, so the small
            % entries disagree in the 17th significant figure -- a relative error of 0.17 on
            % a number that is zero for every purpose. Scale the tolerance to the ARRAY.
            %
            % ⚠ AND THE REFERENCE MUST BE THE UNCANCELLED MAGNITUDE. S1/S2 are pooled moduli,
            % so their sums cannot cancel and the array max is the right reference. S0 is a
            % SIGNED mean and cancels completely (see s0VanishesOnBandpass), so scaling to
            % max|S0| would be scaling to zero. Pass sum|w||F| for that one -- the backward
            % error of a dot product goes as sum of |terms|, not as |sum of terms|.
            t = k * max(abs(ref(:)));
        end

        function F = field(fx, nT, seed)
            % A complex field, as the pipeline actually supplies it: a few localized atoms
            % with independent complex time courses. Complex on purpose -- |.| is then a
            % genuine envelope rather than a rectification (see the scatter docstring).
            rng(seed);
            seeds = [7 700 1500];
            A = rheome.filters.localize(fx.basis, seeds(:), fx.frame.g{2}(fx.basis.Lambda));
            A = A ./ max(abs(A), [], 1);                 % unit-peak atoms
            F = A * (randn(numel(seeds), nT) + 1i*randn(numel(seeds), nT));
        end
    end

    methods (Test)

        function shapes(tc)
            fx = tc.fixture();  nT = 12;
            S  = rheome.filters.scatter(fx.basis, tc.field(fx, nT, 1), fx.frame, ...
                                 struct('MaxBytes', tc.MB));
            M = fx.M;
            tc.verifySize(S.S0, [1 nT]);
            tc.verifySize(S.S1, [M 1 nT]);
            tc.verifySize(S.S2, [M M 1 nT]);
            tc.verifyEqual(S.Order, 2);
            tc.verifyTrue(all(isfinite(S.S1(:))), 'S1 must be finite everywhere.');
        end

        function nonnegative(tc)
            % Modulus, then non-negative pooling weights. Nothing can make these negative.
            fx = tc.fixture();
            S  = rheome.filters.scatter(fx.basis, tc.field(fx, 8, 2), fx.frame, ...
                                 struct('MaxBytes', tc.MB));
            tc.verifyGreaterThanOrEqual(S.S1(:), 0);
            v = S.S2(~isnan(S.S2));
            tc.verifyGreaterThanOrEqual(v, 0);
            tc.verifyNotEmpty(v, 'the second order computed no pairs at all.');
        end

        function homogeneous(tc)
            % Scattering is 1-homogeneous: every stage is linear except |.|, and |a z| =
            % a|z| for real a > 0. This is exact, so it is asserted at machine precision.
            fx = tc.fixture();  F = tc.field(fx, 6, 3);  a = 3.7;
            o  = struct('MaxBytes', tc.MB);
            S1 = rheome.filters.scatter(fx.basis, F,     fx.frame, o);
            S2 = rheome.filters.scatter(fx.basis, a * F, fx.frame, o);
            tc.verifyEqual(S2.S1, a * S1.S1, 'AbsTol', tc.tol(a * S1.S1, 1e-12));
            k  = ~isnan(S1.S2);
            tc.verifyEqual(S2.S2(k), a * S1.S2(k), 'AbsTol', tc.tol(a * S1.S2(k), 1e-12));
        end

        function blockInvariant(tc)
            % ⭐ THE MEMORY BOUND MUST NOT CHANGE THE ANSWER. Time blocking is the whole
            % reason this runs at cortex scale; if a small budget gave a different result
            % the guard would be silently corrupting every large run.
            fx = tc.fixture();  F = tc.field(fx, 17, 4);
            whole = rheome.filters.scatter(fx.basis, F, fx.frame, struct('MaxBytes', tc.MB));
            tiny  = rheome.filters.scatter(fx.basis, F, fx.frame, struct('MaxBytes', 1));
            tc.verifyEqual(tiny.BlockSize, 1);
            tc.verifyEqual(tiny.NumBlocks, 17);
            tc.verifyGreaterThan(whole.BlockSize, 1);
            a = full(sum(fx.basis.Mass, 2));  a = a / sum(a);
            tc.verifyEqual(tiny.S0, whole.S0, 'AbsTol', tc.tol(a.' * abs(F), 1e-12));
            tc.verifyEqual(tiny.S1, whole.S1, 'AbsTol', tc.tol(whole.S1, 1e-12));
            k = ~isnan(whole.S2);
            tc.verifyEqual(isnan(tiny.S2), isnan(whole.S2));
            tc.verifyEqual(tiny.S2(k), whole.S2(k), 'AbsTol', tc.tol(whole.S2(k), 1e-12));
        end

        function coarserSecondRule(tc)
            % Mallat's lambda2 < lambda1, read on sigma: a pair is computed iff the SECOND
            % member is strictly coarser. Everything else stays NaN -- not zero, because
            % "not computed" and "computed and came out zero" must stay distinguishable.
            fx = tc.fixture();
            S  = rheome.filters.scatter(fx.basis, tc.field(fx, 5, 5), fx.frame, ...
                                 struct('MaxBytes', tc.MB));
            sig = S.Sigma;  M = fx.M;
            for m1 = 1:M
                for m2 = 1:M
                    got = squeeze(S.S2(m1, m2, 1, :));
                    if sig(m2) > sig(m1)
                        tc.verifyTrue(all(~isnan(got)), ...
                            sprintf('pair (%d,%d) should be computed.', m1, m2));
                    else
                        tc.verifyTrue(all(isnan(got)), ...
                            sprintf('pair (%d,%d) should be excluded.', m1, m2));
                    end
                end
            end
            tc.verifyEqual(size(S.Pairs,1), sum(sum(sig(:).' > sig(:))));
        end

        function constantFieldFallsThrough(tc)
            % A constant field is the lambda = 0 mode alone. Every band-pass member has
            % g(0) = 0, so S1 must vanish on them -- and NOT on the low-pass member, which
            % is what tells the two kinds of member apart.
            fx = tc.fixture();
            F  = ones(fx.nV, 3);
            S  = rheome.filters.scatter(fx.basis, F, fx.frame, ...
                                 struct('MaxBytes', tc.MB, 'Order', 1));
            band = isfinite(S.Sigma);
            tc.verifyLessThan(max(reshape(S.S1(band,:,:), [], 1)), 1e-9);
            tc.verifyGreaterThan(max(reshape(S.S1(~band,:,:), [], 1)), 1e-3);
        end

        function matchedScalePeaks(tc)
            % Plant an atom built from member m and check S1 peaks on m -- the reason a
            % member index can be read as a physical size at all.
            %
            % ⚠ THE BAND IS BOUNDED AT BOTH ENDS, AND ONLY THE MIDDLE IS ASSERTED. Measured
            % on this fixture (sigma in mm, atom -> S1 argmax, identical for seeds 700/1500):
            %
            %   m        2     3     4     5  |  6     7     8     9  |  10    11    12    13
            %   sigma  11.0  13.9  17.5  22.1 | 27.9  35.1  44.3  55.9| 70.5  88.9 112.1 141.3
            %   peak      2     2/3   3     4 |  6     7     8     9  |   9    11    11    11
            %
            % Below m = 6 the members sit under frame.SigmaFloor (24.0 mm): the basis cannot
            % resolve them, they under-respond, and the peak slides. Above m = 9 sigma passes
            % ~0.6 R on a 100 mm sphere, the atom wraps the domain and coarse members saturate
            % onto each other. Neither end is a defect in scatter -- both are the surface and
            % the basis, and asserting across them would be asserting an artefact.
            %
            % ⭐ AND IT IS THE POOLED MODULUS THAT IS MATCHED, NOT THE ENERGY. Over the same
            % band, argmax of the L2 energy sum|g_m c|^2 lands one member FINER than the atom
            % (5->4, 6->5, 7->6), while S1 = <|g_m F|> lands exactly on it. So S1's member
            % index is the trustworthy scale readout, and a spectral-energy proxy is not.
            fx = tc.fixture();
            usable = 6:9;
            tc.assumeGreaterThanOrEqual(numel(fx.frame.g), max(usable));
            for m = usable
                A = rheome.filters.localize(fx.basis, 700, fx.frame.g{m}(fx.basis.Lambda));
                S = rheome.filters.scatter(fx.basis, complex(A), fx.frame, ...
                                    struct('MaxBytes', tc.MB, 'Order', 1));
                [~, peak] = max(S.S1(:, 1, 1));
                tc.verifyEqual(peak, m, sprintf( ...
                    'atom from member %d (%.1f mm) peaked on member %d (%.1f mm).', ...
                    m, 1000*fx.frame.Sigma(m), peak, 1000*fx.frame.Sigma(peak)));
            end
        end

        function s0VanishesOnBandpass(tc)
            % ⭐ THE POOLED MEAN *IS* THE lambda = 0 PROJECTION. <F> = w'F with w the vertex
            % areas is, up to normalisation, F's constant-mode coefficient. A band-pass member
            % has g(0) = 0, so anything it generates has NO constant mode and S0 is exactly
            % zero -- measured 4.9e-17 against field values of 2.6, a cancellation of 1/eps.
            % Worth pinning: it is why S0 carries no information for a band-limited field, and
            % why comparing S0 to itself needs an uncancelled tolerance reference.
            fx = tc.fixture();  F = tc.field(fx, 6, 8);
            S  = rheome.filters.scatter(fx.basis, F, fx.frame, ...
                                 struct('MaxBytes', tc.MB, 'Order', 1));
            a  = full(sum(fx.basis.Mass, 2));  a = a / sum(a);
            tc.verifyLessThan(max(abs(S.S0)), 1e-12 * max(a.' * abs(F)));
        end

        function orderOneSkipsSecond(tc)
            fx = tc.fixture();
            S  = rheome.filters.scatter(fx.basis, tc.field(fx, 4, 6), fx.frame, ...
                                 struct('MaxBytes', tc.MB, 'Order', 1));
            tc.verifyEmpty(S.S2);
            tc.verifyEqual(S.Order, 1);
        end

        function poolingWeightsHonoured(tc)
            % Two disjoint regions must pool independently, and the mass-weighted global
            % mean must be their area-weighted combination.
            fx = tc.fixture();  nV = fx.nV;
            a  = full(sum(fx.basis.Mass, 2));
            half = false(nV,1);  half(1:round(nV/2)) = true;
            P = [half .* a, ~half .* a];
            P = P ./ sum(P, 1);                       % two region means
            F = tc.field(fx, 4, 7);
            o = struct('MaxBytes', tc.MB, 'Order', 1);
            Sg = rheome.filters.scatter(fx.basis, F, fx.frame, o);
            Sr = rheome.filters.scatter(fx.basis, F, fx.frame, setfield(o, 'Pool', P)); %#ok<SFLD>
            tc.verifySize(Sr.S1, [fx.M 2 4]);
            wA = [sum(a(half)); sum(a(~half))] / sum(a);
            recombined = squeeze(sum(Sr.S1 .* reshape(wA, 1, 2, 1), 2));
            tc.verifyEqual(recombined, squeeze(Sg.S1), ...
                'AbsTol', tc.tol(squeeze(Sg.S1), 1e-10));
        end

    end
end

% Author: Diellor Basha, 2026
