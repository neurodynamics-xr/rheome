classdef tSensitivity < matlab.unittest.TestCase
% rheome.flow.sensitivity -- what a flow kernel produces from NOISE, per vertex.
%
% The normaliser that makes regions comparable. A measured energy is not interpretable on its
% own because the kernel's sensitivity varies enormously over the surface: measured on the
% AnphySleep EEG kernel, the top 1% of vertices carry 15% of the total squared row norm, and
% dividing ROI energy by AREA then inflates small scouts further -- which produced
% 'transversetemporal R' as the peak ROI for every band in every sleep stage.
%
% Author: Diellor Basha, 2026

    properties
        V = 40; C = 6;
        K                       % a vertex operator with deliberately uneven sensitivity
        Cn                      % a noise covariance
    end

    methods (TestClassSetup)
        function build(tc)
            rng(4);
            tc.K = randn(tc.V, tc.C);
            tc.K(1:5, :) = tc.K(1:5, :) * 10;      % five hot vertices
            A = randn(tc.C);  tc.Cn = A*A' + eye(tc.C);
        end
    end

    methods (Test)

        function equalsTheVarianceOfTheKernelOutputUnderThatNoise(tc)
            % s_v = k_v * Cn * k_v' -- the variance the kernel produces at v from noise.
            s = rheome.flow.sensitivity(tc.K, tc.Cn);
            tc.verifySize(s, [tc.V 1]);
            for v = [1 7 tc.V]
                tc.verifyEqual(s(v), tc.K(v,:) * tc.Cn * tc.K(v,:).', 'RelTol', 1e-12);
            end
        end

        function isNonNegative(tc)
            % A variance cannot be negative; a covariance that is only PSD to rounding can
            % make it so, and a negative normaliser would flip the sign of every ratio.
            s = rheome.flow.sensitivity(tc.K, tc.Cn);
            tc.verifyTrue(all(s >= 0));
        end

        function identityCovarianceGivesTheSquaredRowNorm(tc)
            s = rheome.flow.sensitivity(tc.K, eye(tc.C));
            tc.verifyEqual(s, sum(tc.K.^2, 2), 'RelTol', 1e-12);
        end

        function itTracksTheHotVerticesItIsMeantToCorrect(tc)
            s = rheome.flow.sensitivity(tc.K, tc.Cn);
            tc.verifyGreaterThan(median(s(1:5)) / median(s(6:end)), 20);
        end

        function normalisingRemovesAPureSensitivityDifference(tc)
            % THE POINT. Two regions carrying the SAME underlying activity but seen through
            % different kernel sensitivity must normalise to the same value -- otherwise the
            % map reports the kernel, not the brain.
            s = rheome.flow.sensitivity(tc.K, tc.Cn);
            b = randn(tc.C, 500);                       % identical drive everywhere
            E = sum((tc.K*b).^2, 2) / size(b,2);        % measured energy per vertex
            r = E ./ s;
            tc.verifyLessThan(std(r(1:5))/mean(r(1:5)),  0.5);
            tc.verifyEqual(mean(r(1:5)), mean(r(6:end)), 'RelTol', 0.35);
            % un-normalised, the hot vertices are wildly larger
            tc.verifyGreaterThan(mean(E(1:5))/mean(E(6:end)), 20);
        end

        function rejectsAMismatchedCovariance(tc)
            tc.verifyError(@() rheome.flow.sensitivity(tc.K, eye(tc.C+1)), 'flow:sensitivity:size');
        end

    end
end

% Author: Diellor Basha, 2026
