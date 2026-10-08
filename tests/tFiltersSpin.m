classdef tFiltersSpin < matlab.unittest.TestCase
% rheome.filters.spin: conjugation by a quaternion FIELD, which is what rheome.filters.steer is not.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function aRotationPreservesLength(t)
            rng(3);  V = randn(50,3);  n = randn(50,3);
            Vq = rheome.filters.spin(V, n, 0.7);
            t.verifyEqual(vecnorm(Vq,2,2), vecnorm(V,2,2), 'RelTol', 1e-12);
        end
        function rotatingAboutTheNormalLeavesTheNormalPartExactlyAlone(t)
            % ⭐⭐ the property that makes this usable on a physiologically MIXED field
            rng(4);  n = randn(80,3);  n = n./vecnorm(n,2,2);
            V = randn(80,3);
            Vq = rheome.filters.spin(V, n, pi/3);
            t.verifyEqual(sum(Vq.*n,2), sum(V.*n,2), 'AbsTol', 1e-12);
        end
        function theTangentialPartTurnsByExactlyTheta(t)
            rng(5);  n = randn(60,3);  n = n./vecnorm(n,2,2);
            V = randn(60,3);  Vt = V - sum(V.*n,2).*n;
            for th = [pi/6 pi/2 2.1]
                Wt = rheome.filters.spin(Vt, n, th);
                a = acos(max(-1,min(1, sum(Vt.*Wt,2)./max(vecnorm(Vt,2,2).*vecnorm(Wt,2,2),eps))));
                t.verifyEqual(a, repmat(th, 60, 1), 'AbsTol', 1e-9, sprintf('theta = %g', th));
            end
        end
        function thetaIsTheFullAngleNotTheHalf(t)
            % ⚠ the half-angle belongs to the quaternion and is applied internally
            n = [0 0 1];  V = [1 0 0];
            Vq = rheome.filters.spin(V, n, pi/2);
            t.verifyEqual(Vq, [0 1 0], 'AbsTol', 1e-12);
        end
        function itAcceptsBothLayouts(t)
            rng(6);  V = randn(20,3);  n = randn(20,3);
            a = rheome.filters.spin(V, n, 0.4);
            b = rheome.filters.spin(reshape(V',[],1), n, 0.4);
            t.verifyEqual(b, reshape(a',[],1), 'RelTol', 1e-12);
        end
        function itDiffersFromSteerWhichIsGlobal(t)
            % ⚠ steer is right multiplication by ONE q0; this is conjugation by a FIELD
            rng(7);  nV = 40;  V = randn(nV,3);
            n = randn(nV,3);  n = n./vecnorm(n,2,2);
            th = linspace(0, pi, nV)';                      % a varying angle: steer cannot do this
            Vq = rheome.filters.spin(V, n, th);
            ang = acos(max(-1,min(1, sum(V.*Vq,2)./max(vecnorm(V,2,2).*vecnorm(Vq,2,2),eps))));
            t.verifyGreaterThan(std(ang), 0.3, 'the rotation must vary with position');
        end
    end
end
