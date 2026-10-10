classdef tVortexAtom < matlab.unittest.TestCase
    % rheome.flow.vortexatom: a spinning, non-stationary vortex that stays inside its band.
    %
    % Author: Diellor Basha, 2026

    properties
        B; g; a
    end

    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try
                t.B = rheome.load.bases(rheomeTestSubject());
                t.g = rheome.operators.gauge(t.B.L.S.Vertices, double(t.B.L.S.Faces), Method="diffusion");
            catch
                t.assumeFail('cached bases for test subject are not present');
            end
            t.a = rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false);
        end
    end

    methods (Test)

        function theAtomIsEntirelyInsideItsBand(t)
            % ⭐ the point of using the bank's own member as the time factor
            t.verifyGreaterThan(t.a.check.inbandFrac, 0.99);
            t.verifyEqual(t.a.fc, 8*sqrt(2), 'RelTol', 0.02, 'member is not the alpha centre');
        end

        function everyVertexOscillatesAtTheMemberCentre(t)
            df = t.a.fs / t.a.nT;                      % one DFT bin
            t.verifyEqual(t.a.check.peakHz, t.a.fc, 'AbsTol', 2*df);
        end

        function theSpinIsTheOscillation(t)
            % ⭐⭐ turns = support * fc, NOT support * an independent spin dial bounded by the
            %     tile width. The two-atom construction capped at ~4 turns; this must beat it.
            t.verifyEqual(t.a.check.turnsInSupport, t.a.check.supportSec*t.a.fc, 'RelTol', 0.01);
            t.verifyGreaterThan(t.a.check.turnsInSupport, 10);
            t.verifyEqual(abs(t.a.check.rotationHz), t.a.fc, 'RelTol', 0.05);
        end

        function theEnvelopeIsLocalisedInTime(t)
            % non-stationarity: it does not fill the record, and it is near the tile's support
            t.verifyLessThan(t.a.check.supportSec, 0.5*t.a.nT/t.a.fs);
            % ⚠ rheome.selection.wavelettile says 1.330 s at alpha. That is a DIFFERENT definition
            %   (energy fraction on a different grid), so this is a loose bracket, not equality.
            t.verifyGreaterThan(t.a.check.supportSec, 0.8);
            t.verifyLessThan(t.a.check.supportSec, 2.0);
        end

        function theFieldIsExactlyRankTwo(t)
            J = t.a.materialise();
            t.verifySize(J, [numel(t.a.JA) t.a.nT]);
            for k = [1 round(t.a.nT/3) t.a.nT]
                t.verifyEqual(J(:,k), t.a.JA*t.a.a(k) + t.a.JB*t.a.b(k), 'AbsTol', 0);
            end
        end

        function theQuadraturePartnerIsThePhaseDialTurnedNinety(t)
            % JB must equal the SAME seed with Phase + pi/2, which is what makes the pair free
            s2 = rheome.flow.seedvortex(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false, ...
                                 Phase=t.a.phase + pi/2);
            % both are normalised by their own max, and |1i*z| = |z|, so they must agree
            t.verifyEqual(t.a.JB, s2.J, 'AbsTol', 1e-12);
        end

        function momentNAmSetsTheTotalNotThePeak(t)
            % ⚠ the trap this parameter exists to close: at 140 mm the total is ~652x the peak
            m = rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false, MomentNAm=10);
            t.verifyEqual(m.check.totalMomentNAm, 10, 'RelTol', 1e-9);
            t.verifyLessThan(m.check.peakPerVertexNAm, 1, 'peak should be far below the total');
            t.verifyGreaterThan(m.check.momentPerUnitPeak, 100);
            % scaling must not touch the shape, the band or the spin
            t.verifyEqual(m.check.inbandFrac, t.a.check.inbandFrac, 'RelTol', 1e-9);
            t.verifyEqual(m.check.rotationHz, t.a.check.rotationHz, 'RelTol', 1e-9);
            % ⭐ BOTH columns scale, or the circular motion would turn elliptic
            r = m.check.totalMomentNAm / t.a.check.totalMomentNAm;
            t.verifyEqual(m.JA, t.a.JA*r, 'RelTol', 1e-9);
            t.verifyEqual(m.JB, t.a.JB*r, 'RelTol', 1e-9);
        end

        function totalMomentGrowsWithScaleAtFixedPeak(t)
            % ⚠ why comparing two scales at equal PEAK compares very different sources
            w = zeros(1,3);  lam = [70 140 267];
            for i = 1:3
                q = rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false, ...
                                    WavelengthMM=lam(i));
                w(i) = q.check.momentPerUnitPeak;
            end
            t.verifyTrue(all(diff(w) > 0), 'total moment per unit peak must grow with scale');
            t.verifyGreaterThan(w(3)/w(1), 10, 'measured ratio is ~18x over 70 -> 267 mm');
        end

        function chiralityFlipsOnlyTheSenseOfRotation(t)
            b = rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false, Chirality=-1);
            t.verifyEqual(b.check.rotationHz, -t.a.check.rotationHz, 'RelTol', 1e-6);
            t.verifyEqual(abs(b.psi), abs(t.a.psi), 'AbsTol', 1e-12, 'envelope changed');
            t.verifyEqual(abs(b.zs), abs(t.a.zs), 'AbsTol', 1e-12, 'spatial factor changed');
            t.verifyEqual(b.check.inbandFrac, t.a.check.inbandFrac, 'RelTol', 1e-6);
        end

        function aDifferentBandPicksADifferentMember(t)
            % ⚠ on 8 s, not 4: a 4 s record has no 5.657 Hz member (it returns 6.727) and the
            %   support comes back truncated at 53% of the record, so the doubling fails there
            %   for two separate reasons. Measured at 8 s the ratio is 1.995.
            hi = rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false, Duration=8);
            lo = rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false, ...
                                 Band=[4 8], Duration=8);
            t.verifyLessThan(lo.fc, hi.fc);
            t.verifyNotEqual(lo.member, hi.member);
            t.verifyEqual(lo.check.supportSec/hi.check.supportSec, 2, 'RelTol', 0.02);
            t.verifyGreaterThan(lo.check.inbandFrac, 0.99);
            t.verifyFalse(lo.check.truncated);
        end

        function aTooShortRecordWarnsRatherThanLying(t)
            % ⚠ the trap the doubling test hit: 2.133 s of support inside a 4 s record
            f = @() rheome.flow.vortexatom(rheomeTestSubject(), Bases=t.B, Gauge=t.g, Check=false, ...
                                    Band=[4 8], Duration=4);
            t.verifyWarning(f, 'flow:vortexatom:truncated');
            w = warning('off', 'flow:vortexatom:truncated');
            c = onCleanup(@() warning(w));
            t.verifyTrue(f().check.truncated);
        end

        function theSensorsSeeTheSameBand(t)
            % ⭐ the check the source-side numbers cannot make: push through the leadfield
            try
                st = rheome.load.study(rheomeTestSubject());
                G  = rheome.forward.leadfield(st, GlobalVertices=t.B.L.gv);
            catch
                t.assumeFail('cached study for test subject is not present');
            end
            Bs = (G*t.a.JA)*t.a.a + (G*t.a.JB)*t.a.b;
            P  = abs(fft(Bs, [], 2)).^2;
            ff = (0:t.a.nT-1)*t.a.fs/t.a.nT;  k = 1:floor(t.a.nT/2);
            inb = ff(k) >= t.a.band(1) & ff(k) <= t.a.band(2);
            frac = median(sum(P(:,k(inb)),2) ./ sum(P(:,k),2));
            [~,ip] = max(P(:,k), [], 2);
            t.verifyGreaterThan(frac, 0.99, 'the sensors do not see the source band');
            t.verifyEqual(median(ff(k(ip))), t.a.fc, 'AbsTol', 2*t.a.fs/t.a.nT);
        end

    end
end

% Author: Diellor Basha, 2026
