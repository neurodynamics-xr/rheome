classdef tScaleCleanwelch < matlab.unittest.TestCase
% TSCALECLEANWELCH  rheome.scale.cleanwelch: the Welch reference periodicflow screens tiles against.
%
% The planted case is a real recording's: a 5 s common-mode artefact well inside the recording, which
% rheome.scale.cleanspan masks (not trims) and which used to inflate the PSD every tile is compared to.
% No cached data.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function noMaskIsPwelchBitForBit(tc)
            % the 627 subjects with no interior bad block must not move by one ulp
            rng(1);  fs = 600;  F = randn(20, 60*fs);  nf = 2^nextpow2(8*fs);
            cs = 301:size(F,2) - 150;  mask = false(1, size(F,2));
            mask([1:300, end-149:end]) = true;               % trimmed edges lie outside the span
            [P, f, I] = rheome.scale.cleanwelch(F, fs, cs, mask, nf);
            [P0, f0] = pwelch(F(:, cs).', hann(nf), nf/2, nf, fs);
            tc.verifyEqual(P, P0);  tc.verifyEqual(f, f0);
            tc.verifyEqual(I.nUsed, I.nSeg);
        end

        function segmentsArePwelchsOwn(tc)
            % averaging every segment by hand is the pwelch average: the masked path differs from
            % the unmasked one only in which segments it keeps
            rng(2);  fs = 600;  F = randn(10, 40*fs);  nf = 2^nextpow2(8*fs);  cs = 1:size(F,2);
            mask = false(1, size(F,2));
            [P0, ~, I] = rheome.scale.cleanwelch(F, fs, cs, mask, nf);
            mask(end) = true;                                % touches no segment of pwelch's
            [P1, ~, I1] = rheome.scale.cleanwelch(F, fs, cs, mask, nf);
            tc.verifyEqual(I1.nUsed, I.nSeg);  tc.verifyEqual(P1, P0);
            hop = nf/2;  P = 0;
            for k = 1:I.nSeg
                P = P + pwelch(F(:, (k-1)*hop + (1:nf)).', hann(nf), 0, nf, fs);
            end
            tc.verifyEqual(P / I.nSeg, P0, 'RelTol', 1e-10);
        end

        function interiorArtefactIsLeftOut(tc)
            % there, a masked artefact at 20-25 s used to put every clean tile at power_ratio
            % 0.03-0.07; with it left out of the average the ratio is ~1 again
            rng(3);  fs = 600;  N = 120*fs;  F = randn(40, N);
            t = 20*fs + 1 : 25*fs;  F(:, t) = F(:, t) + 60 * randn(1, numel(t));
            C = rheome.scale.cleanspan(F, fs);
            tc.assertTrue(all(C.mask(t)));
            tc.assertEqual([C.first C.last], [1 N]);          % interior, so masked and not trimmed
            nf = 2^nextpow2(8*fs);  cs = C.first:C.last;
            [P0, f] = pwelch(F(:, cs).', hann(nf), nf/2, nf, fs);
            [P, ~, I] = rheome.scale.cleanwelch(F, fs, cs, C.mask, nf);
            tc.verifyLessThan(I.nUsed, I.nSeg);
            tc.verifyGreaterThan(I.nUsed, 0.75 * I.nSeg);
            nT = 6*fs;  W = fft(F(:, 60*fs + (1:nT)), [], 2);  fw = (0:nT-1)*fs/nT;  pos = 2:nT/2;
            in = fw(pos) >= 1 & fw(pos) <= 45;
            pr = @(Q) mean(abs(W(:, pos(in))).^2, 'all') / ...
                      mean(interp1(f, Q, fw(pos(in))')', 'all') / (fs*nT/2);
            tc.verifyLessThan(pr(P0), 1/3);                   % the old reference rejects a clean tile
            tc.verifyEqual(pr(P), 1, 'AbsTol', 0.15);         % the clean one accepts it
        end

        function allSegmentsMaskedIsAnError(tc)
            rng(4);  fs = 600;  F = randn(5, 20*fs);  mask = true(1, size(F,2));
            tc.verifyError(@() rheome.scale.cleanwelch(F, fs, 1:size(F,2), mask, 4096), 'scale:cleanwelch:nosegments');
        end
    end
end

% Author: Diellor Basha, 2026
