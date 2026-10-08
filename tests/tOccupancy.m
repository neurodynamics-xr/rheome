classdef tOccupancy < matlab.unittest.TestCase
% rheome.detect.onstate (bimodal on/off) and rheome.detect.occupancy (the dyadic occupancy pyramid).
%
% What is asserted: occupancy rolls up exactly; a long episode with short dips is FRAGMENTED at the
% finest level and WHOLE at coarse ones -- the reason for the pyramid; .depth reads the long episode
% as long and a lone short burst as short; and the mixture threshold sits between two planted states
% and reports a single-mode series as not bimodal.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function occupancyRollsUpExactly(tc)
            rng(1);  b = rand(3, 1024) > 0.6;
            O = rheome.detect.occupancy(b, 0.1);
            for L = 1:numel(O.occ) - 1
                c = O.occ{L};  p = O.occ{L+1};  m = size(p, 2);
                tc.verifyEqual(p, (c(:, 1:2:2*m) + c(:, 2:2:2*m)) / 2, 'AbsTol', 1e-15);
            end
            tc.verifyEqual(O.occ{end}, mean(double(b(:, 1:2^(numel(O.occ)-1))), 2), 'AbsTol', 1e-15);
        end

        function aLongEpisodeWithDipsIsWholeAtCoarseLevels(tc)
            dt = 0.1;  n = 1024;  b = false(1, n);
            ep = 101:400;  b(ep) = true;                         % 30 s episode
            for s = ep(1):20:ep(end), b(s:s+1) = false; end      % a 0.2 s dip every 2 s
            b(700:704) = true;                                   % and a lone 0.5 s burst
            O = rheome.detect.occupancy(b, dt, Theta=0.75);
            tc.verifyLessThan(max(O.runs{1}{1}), 2.01);          % finest level: fragments of <= 2 s
            L4 = find(O.tileS >= 1.6, 1);                        % a level with 1.6 s tiles
            tc.verifyGreaterThan(max(O.runs{L4}{1}), 25);        % ... sees the ~30 s episode whole
            tc.verifyGreaterThanOrEqual(median(O.depth(ep(b(ep)))), 12.8);
            tc.verifyLessThanOrEqual(max(O.depth(700:704)), 0.8);
        end

        function mixtureThresholdSitsBetweenPlantedStates(tc)
            rng(2);  y = [exp(randn(1, 3000) * 0.3) exp(2 + randn(1, 2000) * 0.3)];
            S = rheome.detect.onstate(y);
            tc.verifyTrue(S.bimodal);
            tc.verifyGreaterThan(S.threshold, exp(0.5));  tc.verifyLessThan(S.threshold, exp(1.5));
            tc.verifyEqual(S.fracOn, 0.4, 'AbsTol', 0.03);
        end

        function aSingleModeIsNotBimodal(tc)
            rng(3);  S = rheome.detect.onstate(exp(randn(1, 5000) * 0.5));
            tc.verifyFalse(S.bimodal);
        end
    end
end

% Author: Diellor Basha, 2026
