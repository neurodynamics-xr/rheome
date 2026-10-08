classdef tMeasurable < matlab.unittest.TestCase
% rheome.forward.measurable -- what fraction of a hypothesised source pattern the array can see.
%
% Ground truth is linear algebra, so most of this is exact: a pattern built from the retained
% right singular vectors is fully measurable, one from the orthogonal complement is not
% measurable at all, and the fraction is a ratio of energies so it cannot leave [0,1] or
% depend on amplitude. The last test pins the TRAP rather than a capability -- a minimum-norm
% reconstruction scores 1 by construction, and that number must never be read as quality.
%
% Author: Diellor Basha, 2026

    methods (Static)
        function fx = fixture()
            % Synthetic gain, deliberately: this function is linear algebra on a row space, and
            % a real leadfield would add anatomy without testing anything extra. Rank is set
            % well below nV so there IS a null space to find.
            persistent c
            if ~isempty(c), fx = c; return; end
            rng(20260826);
            nV = 300; nCh = 40;
            G  = randn(nCh, 3*nV);
            O  = randn(nV, 3);  O = O ./ vecnorm(O, 2, 2);
            hm = struct('Gain', G, 'GridOrient', O, 'nV', nV, 'nCh', nCh);
            fx = struct('hm', hm, 'nV', nV, 'nCh', nCh, ...
                        'area', 0.5 + rand(nV,1));      % non-uniform, so the metric is exercised
            c = fx;
        end

        function [Vr, sw, r] = rowspace(fx, cut)
            hm = fx.hm;  nV = fx.nV;
            O  = hm.GridOrient ./ vecnorm(hm.GridOrient, 2, 2);
            Osp = sparse((1:3*nV).', repelem((1:nV).',3,1), reshape(O.',[],1), 3*nV, nV);
            sw = sqrt(fx.area);
            Ao = (hm.Gain * Osp) .* sw.';
            [~, sv, Vs] = svd(Ao, 'econ');  sv = diag(sv);
            r  = max(1, sum(sv > sv(1)*cut));
            Vr = Vs(:, 1:r);
        end
    end

    methods (Test)

        function inRowSpaceIsFullyMeasurable(tc)
            fx = tc.fixture();  cut = 1e-3;
            [Vr, sw] = tc.rowspace(fx, cut);
            X = (Vr * randn(size(Vr,2), 4)) ./ sw;          % built from the row space
            m = rheome.forward.measurable(fx.hm, X, struct('Area', fx.area, 'SVCut', cut));
            tc.verifyEqual(m.fraction, ones(4,1), 'AbsTol', 1e-10);
        end

        function inNullSpaceIsInvisible(tc)
            % ⭐ THE POINT OF THE WHOLE FUNCTION. These patterns produce NO field at all: add
            % any of them to a source and not one sensor sample changes.
            fx = tc.fixture();  cut = 1e-3;
            [Vr, sw, r] = tc.rowspace(fx, cut);
            Z = randn(fx.nV, 4);
            Z = Z - Vr * (Vr.' * Z);                        % orthogonal complement
            X = Z ./ sw;
            m = rheome.forward.measurable(fx.hm, X, struct('Area', fx.area, 'SVCut', cut));
            tc.verifyLessThan(max(m.fraction), 1e-10);
            tc.verifyEqual(m.rank, r);
            tc.verifyLessThan(r, fx.nV, 'the fixture must have a null space to find.');
        end

        function fractionIsBoundedAndAmplitudeFree(tc)
            fx = tc.fixture();
            X  = randn(fx.nV, 6);
            o  = struct('Area', fx.area);
            m1 = rheome.forward.measurable(fx.hm, X,       o);
            m2 = rheome.forward.measurable(fx.hm, 137.5*X, o);
            tc.verifyGreaterThanOrEqual(m1.fraction, 0);
            tc.verifyLessThanOrEqual(m1.fraction, 1 + 1e-12);
            tc.verifyEqual(m2.fraction, m1.fraction, 'RelTol', 1e-12);
        end

        function projectionIsIdempotent(tc)
            fx = tc.fixture();  o = struct('Area', fx.area, 'Project', true);
            m1 = rheome.forward.measurable(fx.hm, randn(fx.nV, 3), o);
            m2 = rheome.forward.measurable(fx.hm, m1.projected,    o);
            tc.verifyEqual(m2.fraction, ones(3,1), 'AbsTol', 1e-10);
        end

        function rankMovesWithTheCut(tc)
            % The fraction is "measurable at this assumed SNR", not an absolute -- a looser cut
            % keeps more weakly-observed directions and can only raise the fraction.
            fx = tc.fixture();  X = randn(fx.nV, 3);
            loose = rheome.forward.measurable(fx.hm, X, struct('Area',fx.area,'SVCut',1e-6));
            tight = rheome.forward.measurable(fx.hm, X, struct('Area',fx.area,'SVCut',1e-1));
            tc.verifyGreaterThanOrEqual(loose.rank, tight.rank);
            tc.verifyGreaterThanOrEqual(loose.fraction + 1e-12, tight.fraction);
        end

        function unconstrainedAccepted(tc)
            fx = tc.fixture();
            m = rheome.forward.measurable(fx.hm, randn(3*fx.nV, 2), struct('Area', fx.area));
            tc.verifySize(m.fraction, [2 1]);
            tc.verifyLessThanOrEqual(m.fraction, 1 + 1e-12);
            tc.verifyLessThanOrEqual(m.rank, fx.nCh);
        end

        function minimumNormScoresOneTautologically(tc)
            % ⚠ PINNING THE TRAP, NOT A CAPABILITY. A min-norm estimate is the solution whose
            % null-space component is exactly zero, so it is fully "measurable" whatever the
            % data was -- including data that is pure noise, as here. Never read this number
            % as evidence that a reconstruction is good.
            fx = tc.fixture();  nV = fx.nV;
            O   = fx.hm.GridOrient ./ vecnorm(fx.hm.GridOrient,2,2);
            Osp = sparse((1:3*nV).', repelem((1:nV).',3,1), reshape(O.',[],1), 3*nV, nV);
            Gc  = fx.hm.Gain * Osp;
            b   = randn(fx.nCh, 5);                                  % pure noise "data"
            xhat = Gc.' * ((Gc*Gc.' + 1e-9*eye(fx.nCh)) \ b);        % minimum norm
            m = rheome.forward.measurable(fx.hm, xhat, struct('Area', ones(nV,1), 'SVCut', 1e-8));
            tc.verifyEqual(m.fraction, ones(5,1), 'AbsTol', 1e-6);
        end

    end
end

% Author: Diellor Basha, 2026
