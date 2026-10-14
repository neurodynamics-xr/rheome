classdef tScaleRamp < matlab.unittest.TestCase
% MS1 2.4: a peak whose size ramps is read by rheome.detect.ridges at its planted size (gamma = 1 for a scalar
% peak), on one ridge, at the planted vertex -- on the in-job sphere and on a non-spherical "cortex".
%
% Author: Diellor Basha, 2026

    methods (Test)
        function readsThePlantedSize(tc)
            [V, F] = rheome.geom.icosphere(4);  V = 0.09 * V .* [1 1 0.8];           % an ellipsoid as the cortex
            [L, M] = rheome.operators.laplace_beltrami(V, F);  b = rheome.eigen.modes(L, M, 800);
            S.B.L = struct('S', struct('Vertices', V, 'Faces', F), 'lbo', b);
            [T, X] = rheome.scale.measure_scaleramp('synth', S, SigmaMM=[10 30], NFrames=33, BankMM=[5 60], ...
                                                    SphereSub=4, KSphere=800);
            v = @(m, s) T.value(T.metric == m & T.band == s);
            for s = ["sphere" "cortex"]
                tc.verifyEqual(v("slope", s), 1, 'AbsTol', 0.1, s);
                tc.verifyEqual(v("ratio_median", s), 1, 'AbsTol', 0.1, s);
                tc.verifyGreaterThanOrEqual(v("continuity", s), 0.9, s);
                tc.verifyLessThan(v("centre_err_mm_max", s), 10, s);
            end
            tc.verifyEqual(height(X), 2 * 33);
        end
    end
end

% Author: Diellor Basha, 2026
