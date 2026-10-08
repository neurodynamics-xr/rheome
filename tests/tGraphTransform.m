classdef tGraphTransform < matlab.unittest.TestCase

    methods (Test)

        function eigenRoundTripsAnInBasisField(tc)
            fx = gfbFixture();
            T  = rheome.graphtransform.eigen(fx.Phi, fx.Mass, fx.Lambda);
            F  = fx.Phi * randn(numel(fx.Lambda), 3);      % lives in the basis exactly
            tc.verifyEqual(T.inverse(T.forward(F)), F, 'AbsTol', 1e-8);
        end

        function eigenForwardIsPhiTransposeMassF(tc)
            fx = gfbFixture();
            T  = rheome.graphtransform.eigen(fx.Phi, fx.Mass, fx.Lambda);
            F  = randn(fx.nV, 2);
            tc.verifyEqual(T.forward(F), fx.Phi' * (fx.Mass * F), 'AbsTol', 1e-12);
        end

        function scalarNormIsAbs(tc)
            fx = gfbFixture();
            T  = rheome.graphtransform.eigen(fx.Phi, fx.Mass, fx.Lambda);
            F  = randn(fx.nV, 4);
            tc.verifyEqual(T.norm(F), abs(F), 'AbsTol', 1e-12);
        end

        function quaternionNormUsesImaginaryRowsOnly(tc)
            % The quaternion special case is a NORM on the field type, and that is the
            % only thing Dirac changes.
            nV  = 5;  K = 4;
            Phi = randn(4*nV, K);  B = speye(4*nV);
            T   = rheome.graphtransform.eigen(Phi, B, (1:K)', 'nV', nV);
            F   = zeros(4*nV, 1);
            F(1:4:end) = 99;                   % real part -- must be ignored
            F(2:4:end) = 3;  F(3:4:end) = 4;   % (x,y) -> magnitude 5
            tc.verifyEqual(T.norm(F), 5*ones(nV,1), 'AbsTol', 1e-12);
        end

        function validateAcceptsAGoodTransform(tc)
            fx = gfbFixture();
            tc.verifyTrue(rheome.graphtransform.validate( ...
                rheome.graphtransform.eigen(fx.Phi, fx.Mass, fx.Lambda)));
        end

        function validateRejectsAMissingHandle(tc)
            [ok, why] = rheome.graphtransform.validate(struct('forward', @(F) F));
            tc.verifyFalse(ok);
            tc.verifyNotEmpty(why);
        end

        function settingABadTransformErrors(tc)
            tc.verifyError(@() rheome.graphfilterbank(100, 'Transform', struct('forward', @(F) F)), ...
                'graphfilterbank:transform');
        end

        function goodTransformAttaches(tc)
            fx = gfbFixture();
            g  = rheome.graphfilterbank(fx.Lambda, ...
                     'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, fx.Lambda));
            tc.verifyNotEmpty(g.Transform);
            tc.verifyEqual(g.Transform.rows, fx.nV);
        end

    end
end

% Author: Diellor Basha, 2026
