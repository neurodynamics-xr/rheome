classdef tJfbImpulse < matlab.unittest.TestCase

    methods (Test)

        function speedKernelPeaksOnItsDiagonal(tc)
            % The kernel selects omega = c*sqrt(lambda).
            c = 5;  K = rheome.jointfilterbank.speedkernel(c, 2);
            lam = (1:10).';  om = c*sqrt(lam).';
            G = K(lam, om);
            [~, ipk] = max(G, [], 2);
            tc.verifyEqual(ipk, (1:10).');       % the peak is on the diagonal
        end

        function speedKernelBroadcasts(tc)
            K = rheome.jointfilterbank.speedkernel(5, 2);
            tc.verifySize(K((1:7).', (1:4)), [7 4]);
        end

        function speedKernelMakesABankNonSeparable(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                     'JointKernels', {rheome.jointfilterbank.speedkernel(5, 20)});
            tc.verifyFalse(j.Separable);
        end

        function impulseReturnsAFieldOverTime(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Half','full', 'Transform', fx.T);
            U = impulse(j, 1, 100);
            tc.verifySize(U, [fx.nV, 64]);
            tc.verifyTrue(isreal(U));
        end

        function impulseIsLocalisedAtItsSeedAtTimeZero(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Half','full', 'Transform', fx.T);
            seed = 100;
            U = impulse(j, 1, seed);
            [~, ipk] = max(abs(U(:,1)));
            d = norm(fx.V(ipk,:) - fx.V(seed,:));
            tc.verifyLessThan(d, 0.05);          % within 50 mm on a 100 mm sphere
        end

        function transformIsInheritedFromTheGraphBank(tc)
            % Under composition the field-domain route belongs to the graph side.
            fx  = jfbFixture();
            gfb = rheome.graphfilterbank(fx.Lambda, 'Wavelet','itersine', 'NumFilters',3, ...
                                  'Transform', fx.T);
            j = rheome.jointfilterbank(gfb, 'SignalLength',64, 'Half','full');
            tc.verifyNotEmpty(j.Transform);
            tc.verifySize(impulse(j, 1, 50), [fx.nV, 64]);
        end

        function anExplicitTransformOverridesTheInherited(tc)
            fx  = jfbFixture();
            gfb = rheome.graphfilterbank(fx.Lambda, 'Wavelet','itersine', 'NumFilters',3, ...
                                  'Transform', fx.T);
            j = rheome.jointfilterbank(gfb, 'SignalLength',64, 'Transform', fx.T);
            tc.verifyNotEmpty(j.Transform);
        end

        function impulseNeedsATransform(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Half','full');
            tc.verifyError(@() impulse(j, 1, 1), 'jointfilterbank:noTransform');
        end

        function impulseNeedsAFullDftGrid(tc)
            % With arbitrary Frequencies there is no bin mapping back to time.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'Frequencies', linspace(8,13,17), ...
                                 'Transform', fx.T);
            tc.verifyError(@() impulse(j, 1, 1), 'jointfilterbank:noGrid');
        end

    end
end

% Author: Diellor Basha, 2026
