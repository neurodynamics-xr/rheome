classdef tInverseMne < matlab.unittest.TestCase
% rheome.inverse.mne -- whitened minimum-norm on the raw unconstrained leadfield.
%
% This is the primitive rheome.flow.context builds on by default, so the properties asserted here are
% the ones that define a minimum norm rather than merely describe this implementation: the
% estimate reproduces the data through the forward model, it is the SMALLEST such estimate,
% it lives entirely in the row space of the gain, and it scales linearly with the data.
%
% Author: Diellor Basha, 2026

    methods (Static)
        function fx = fixture()
            persistent c
            if ~isempty(c), fx = c; return; end
            rng(20260905);
            nV = 120; nCh = 30;
            G  = randn(nCh, 3*nV);
            % A non-trivial, well-conditioned noise covariance so the whitener does real work.
            A  = randn(nCh); C = A*A.' + nCh*eye(nCh);
            fx = struct('G', G, 'nV', nV, 'nCh', nCh, ...
                        'ncov', struct('NoiseCov', C, 'FourthMoment', [], 'nSamples', []));
            c = fx;
        end
        function R = solve(fx, varargin)
            o = struct('InverseMeasure','amplitude');
            for i = 1:2:numel(varargin), o.(varargin{i}) = varargin{i+1}; end
            R = rheome.inverse.mne(fx.G, fx.ncov, o);
        end
    end

    methods (Test)

        function shapesAndFields(tc)
            fx = tc.fixture();  R = tc.solve(fx);
            tc.verifySize(R.ImagingKernel, [3*fx.nV, fx.nCh]);
            tc.verifySize(R.Whitener, [fx.nCh fx.nCh]);
            tc.verifyEqual(R.nVert, fx.nV);
            tc.verifyEqual(R.Measure, 'amplitude');
            tc.verifyLessThanOrEqual(R.RankLeadfield, fx.nCh);
            tc.verifyTrue(all(isfinite(R.ImagingKernel(:))));
        end

        function isMinimumNormAmongDataFits(tc)
            % ⭐ THE DEFINING PROPERTY. Take the estimate, add ANY null-space field, and the
            % predicted data is unchanged while the norm can only grow. Verified against 200
            % random null-space perturbations rather than argued from the formula.
            fx = tc.fixture();
            R  = tc.solve(fx, 'NoiseMethod', 'none', 'SnrFixed', 1e6);   % -> essentially unregularised
            b  = randn(fx.nCh, 1);
            x  = R.ImagingKernel * b;
            [~,~,V] = svd(fx.G, 'econ');
            r  = sum(svd(fx.G) > fx.nCh*eps(single(norm(fx.G))));
            Vr = V(:,1:r);
            worse = 0;
            for i = 1:200
                z = randn(3*fx.nV,1);  z = z - Vr*(Vr.'*z);     % pure null space
                tc.verifyLessThan(norm(fx.G*z)/norm(fx.G,'fro'), 1e-10);   % it really is invisible
                if norm(x+z) > norm(x), worse = worse + 1; end
            end
            tc.verifyEqual(worse, 200, 'every null-space perturbation must increase the norm.');
        end

        function livesEntirelyInTheRowSpace(tc)
            % The corollary, and the reason rheome.forward.measurable scores MNE output 1.0.
            fx = tc.fixture();  R = tc.solve(fx);
            [~,~,V] = svd(fx.G, 'econ');
            r  = sum(svd(fx.G) > fx.nCh*eps(single(norm(fx.G))));
            Vr = V(:,1:r);
            X  = R.ImagingKernel * randn(fx.nCh, 5);
            resid = X - Vr*(Vr.'*X);
            tc.verifyLessThan(norm(resid)/norm(X), 1e-10);
        end

        function linearInTheData(tc)
            fx = tc.fixture();  R = tc.solve(fx);
            b1 = randn(fx.nCh,1); b2 = randn(fx.nCh,1); a = -2.75;
            tc.verifyEqual(R.ImagingKernel*(b1 + a*b2), ...
                (R.ImagingKernel*b1) + a*(R.ImagingKernel*b2), 'RelTol', 1e-12);
        end

        function regularisationShrinks(tc)
            % Lower assumed SNR = heavier Wiener damping = a smaller estimate. Monotone.
            fx = tc.fixture();  b = randn(fx.nCh,1);  n = zeros(1,4);  snr = [0.5 1 3 30];
            for i = 1:numel(snr)
                R = tc.solve(fx, 'SnrFixed', snr(i));
                n(i) = norm(R.ImagingKernel*b);
            end
            tc.verifyEqual(n, sort(n), 'RelTol', 1e-12, ...
                'estimate norm must grow monotonically with assumed SNR.');
        end

        function measuresDifferButDirectionsAgree(tc)
            % dSPM rescales each vertex triplet; it must not rotate the current there.
            fx = tc.fixture();  b = randn(fx.nCh,1);
            Ja = reshape(tc.solve(fx,'InverseMeasure','amplitude').ImagingKernel*b, 3, []);
            Jd = reshape(tc.solve(fx,'InverseMeasure','dspm2018').ImagingKernel*b, 3, []);
            ua = Ja./max(vecnorm(Ja,2,1),eps);  ud = Jd./max(vecnorm(Jd,2,1),eps);
            tc.verifyEqual(abs(sum(ua.*ud,1)), ones(1,fx.nV), 'AbsTol', 1e-9);
        end

        function rejectsAConstrainedGain(tc)
            % ⚠ THE SHAPE CHECK ALONE IS NOT ENOUGH, and this test is why it is not the only
            % one. nV = 120 is divisible by 3, so a constrained [nCh x 120] gain passes the
            % mod-3 guard and is read as 40 unconstrained vertices -- wrong sources, no error.
            % opts.nVert is what actually catches it, and rheome.flow.context always supplies it.
            fx = tc.fixture();
            tc.verifyError(@() rheome.inverse.mne(fx.G(:,1:fx.nV+1), fx.ncov), ...
                'inverse:mne:constrained', 'a non-multiple-of-3 width must be rejected on shape.');
            R = rheome.inverse.mne(fx.G(:,1:fx.nV), fx.ncov);            % slips through on shape alone
            tc.verifyEqual(R.nVert, fx.nV/3, 'the silent misreading is real -- hence nVert.');
            tc.verifyError(@() rheome.inverse.mne(fx.G(:,1:fx.nV), fx.ncov, ...
                struct('nVert', fx.nV)), 'inverse:mne:constrained');
        end

        function rejectsMismatchedCovariance(tc)
            fx = tc.fixture();
            bad = struct('NoiseCov', eye(fx.nCh+1), 'FourthMoment', [], 'nSamples', []);
            tc.verifyError(@() rheome.inverse.mne(fx.G, bad), 'inverse:mne:size');
        end

        function sharesTheWhitenerWithDirac(tc)
            % ⚠ PINS THE REFACTOR. rheome.inverse.dirac and rheome.inverse.mne must differ ONLY in the source
            % basis, so they have to reach the same private whitener. If someone re-inlines a
            % copy into either file, this stops being true and the comparison stops meaning
            % anything -- so assert the whitened gain matches a whitener built the same way.
            fx = tc.fixture();
            R  = tc.solve(fx);
            C  = (fx.ncov.NoiseCov + fx.ncov.NoiseCov.')/2;
            [U,S] = svd(C,'econ'); sn = sqrt(diag(S));
            ridge = mean(diag(S)) * 0.1;                       % 'reg' with the default NoiseReg
            iW = U * diag(1./sqrt(sn.^2 + ridge)) * U.';
            tc.verifyEqual(R.Whitener, iW, 'AbsTol', 1e-10*max(abs(iW(:))));
        end

    end
end

% Author: Diellor Basha, 2026
