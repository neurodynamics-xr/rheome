classdef tPhaseLock < matlab.unittest.TestCase
% rheome.flow.phaselock -- reproducible organisation of a field by a reference rhythm's phase.
% Author: Diellor Basha, 2026

    methods (Static)
        function [z, ph] = lockedField(nV, nCyc, fs, f0, lockFrac, seed)
            rng(seed);
            t  = (0:round(nCyc*fs/f0)-1)/fs;
            ph = angle(exp(1i*2*pi*f0*t));
            pat = randn(nV,1);                                  % a fixed spatial pattern
            zl = pat * exp(1i*2*pi*f0*t);                       % locked to the reference
            zn = (randn(nV,numel(t)) + 1i*randn(nV,numel(t)));  % unlocked
            z  = lockFrac*zl + (1-lockFrac)*zn;
        end
    end

    methods (Test)

        function aPerfectlyLockedFieldReproducesAcrossHalves(tc)
            [z, ph] = tPhaseLock.lockedField(50, 24, 600, 10, 1.0, 1);
            out = rheome.flow.phaselock(z, ph, 30);
            tc.verifyGreaterThan(median(out.splitHalf), 0.95);
            tc.verifyGreaterThan(out.global, 0.95);
            tc.verifyGreaterThan(median(out.plv), 0.9);
        end

        function anUnlockedFieldDoesNotReproduce(tc)
            [z, ph] = tPhaseLock.lockedField(50, 24, 600, 10, 0.0, 2);
            out = rheome.flow.phaselock(z, ph, 30);
            tc.verifyLessThan(abs(median(out.splitHalf)), 0.2);
            tc.verifyLessThan(abs(out.global), 0.2);
        end

        function plvIsBiasedByCycleCountButSplitHalfIsNot(tc)
            % THE REASON splitHalf IS THE HEADLINE. With no locking at all, PLV rises as the
            % cycle count falls -- ~1/sqrt(C) -- and at a handful of cycles it looks like a
            % real effect. The split-half statistic stays at 0 however few cycles there are.
            plv = zeros(1,2);  sh = zeros(1,2);
            nc  = [6 60];
            for k = 1:2
                [z, ph] = tPhaseLock.lockedField(60, nc(k), 600, 10, 0.0, 10+k);
                out = rheome.flow.phaselock(z, ph, 30);
                plv(k) = median(out.plv);  sh(k) = abs(median(out.splitHalf));
            end
            tc.verifyGreaterThan(plv(1), 1.6*plv(2));           % PLV inflates at few cycles
            tc.verifyLessThan(max(sh), 0.25);                   % split-half does not
        end

        function splittingBySampleWouldInflateSoItSplitsByCycle(tc)
            % Neighbouring samples inside one cycle are near-duplicates. The implementation
            % must alternate whole cycles; if it split by sample, an UNLOCKED field would
            % still reproduce strongly because each half holds the other's neighbours.
            [z, ph] = tPhaseLock.lockedField(60, 40, 600, 10, 0.0, 5);
            out = rheome.flow.phaselock(z, ph, 30);
            tc.verifyLessThan(abs(out.global), 0.25);
            tc.verifyGreaterThan(out.nCycle, 30);
        end

        function aPartiallyLockedFieldLandsBetween(tc)
            [z, ph] = tPhaseLock.lockedField(60, 40, 600, 10, 0.35, 7);
            out = rheome.flow.phaselock(z, ph, 30);
            tc.verifyGreaterThan(out.global, 0.3);
            tc.verifyLessThan(out.global, 0.98);
        end

        function standingAndTravellingPatternsAreToldapart(tc)
            % A rank test CANNOT do this: at one frequency the binned map is rank 1 whatever
            % the pattern does. The split has to be read off the spatially-varying PHASE.
            rng(31);
            nV = 200;  fs = 600;  f0 = 10;  t = (0:round(40*fs/f0)-1)/fs;
            ph = angle(exp(1i*2*pi*f0*t));
            pat = randn(nV,1);
            zS = (pat)              * exp(1i*2*pi*f0*t);          % one common phase
            zT = (pat.*exp(2i*pi*(1:nV).'/nV)) * exp(1i*2*pi*f0*t); % phase ramps across space
            oS = rheome.flow.phaselock(zS, ph, 30);
            oT = rheome.flow.phaselock(zT, ph, 30);
            tc.verifyGreaterThan(oS.standing, 0.98);              % standing
            tc.verifyLessThan(   oT.standing, 0.60);              % travelling
            tc.verifyGreaterThan(oS.standing, oT.standing + 0.35);
        end

        function aRealFieldIsRejected(tc)
            [z, ph] = tPhaseLock.lockedField(20, 20, 600, 10, 1.0, 3);
            tc.verifyError(@() rheome.flow.phaselock(real(z), ph, 30), 'flow:phaselock:real');
        end

    end
end
