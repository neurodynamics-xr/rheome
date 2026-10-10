% Author: Diellor Basha, 2026
classdef tQcRules < matlab.unittest.TestCase
% TQCRULES  Detector v3: PSD_HIGH (65–115 Hz) never labels a channel on its own, only corroborates;
% every other rule labels alone, exactly as in v1/v2.
%
% See also: rheome.qc.badchannelrules, rheome.qc.badchannels

methods (Test)
    function eachRuleAloneLabelsExceptPsdHigh(tc)
        [M, P] = tQcRules.i_case();
        [flag, why, rawFlag, rawWhy] = rheome.qc.badchannelrules(M, P);
        tc.verifyEqual(flag', logical([0 1 1 1 1 0 1 1]));
        tc.verifyEqual(rawFlag', logical([0 1 1 1 1 1 1 1]));
        tc.verifyEqual(why([2 3 4 5 6]), {'FLAT'; 'NOISY'; 'JUMPY'; 'PSD_LOW'; ''});
        tc.verifyEqual(rawWhy{6}, 'PSD_HIGH');         % kept for the record, not labelled
    end
    function psdHighCorroborates(tc)
        [M, P] = tQcRules.i_case();
        [~, why] = rheome.qc.badchannelrules(M, P);
        tc.verifyEqual(why{7}, 'JUMPY PSD_HIGH');
        tc.verifyEqual(why{8}, 'PSD_LOW PSD_HIGH');
    end
end

methods (Static)
    function [M, P] = i_case()
        % 8 channels: clean, FLAT, NOISY, JUMPY, PSD_LOW, PSD_HIGH alone, JUMPY+PSD_HIGH, PSD_LOW+PSD_HIGH
        P = struct('zFlat', -5, 'rFlat', 0.2, 'zNoisy', 5, 'rNoisy', 3, 'jumpFrac', 0.2, 'zPsd', 5, 'lrPsd', 0.5);
        o = ones(8, 1);
        M = struct('sd', o, 'z_sd', 0*o, 'r_sd', o, 'jumpfrac', 0*o, 'z_lo', 0*o, 'r_lo', o, 'z_hi', 0*o, 'r_hi', o);
        M.sd(2) = 0;
        M.z_sd(3) = 8; M.r_sd(3) = 5;
        M.jumpfrac([4 7]) = 0.5;
        M.z_lo([5 8]) = 7; M.r_lo([5 8]) = 10;
        M.z_hi([6 7 8]) = 7; M.r_hi([6 7 8]) = 5;
    end
end
end
% Author: Diellor Basha, 2026
