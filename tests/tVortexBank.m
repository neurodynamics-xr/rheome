classdef tVortexBank < matlab.unittest.TestCase
    % rheome.flow.vortexbank / rheome.flow.vortexmatch: what a vortex dictionary can and cannot settle.
    %
    % Author: Diellor Basha, 2026

    properties
        B; g; st; bank; G
    end

    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try
                t.B  = rheome.load.bases(rheomeTestSubject());
                t.g  = rheome.operators.gauge(t.B.L.S.Vertices, double(t.B.L.S.Faces), Method="diffusion");
                t.st = rheome.load.study(rheomeTestSubject());
            catch
                t.assumeFail('cached bases/study for test subject are not present');
            end
            t.G = rheome.forward.leadfield(t.st, GlobalVertices=double(t.B.L.gv(:))');
            t.bank = rheome.flow.vortexbank(rheomeTestSubject(), Locations=16, Bases=t.B, Gauge=t.g, ...
                                     Study=t.st, Verbose=false);
        end
    end

    methods (Access = private)
        function [D, JA, JB] = i_plant(t, vtx, lam, kind)
            sv = rheome.flow.seedvortex(rheomeTestSubject(), Vertex=vtx, WavelengthMM=lam, ...
                                 Bases=t.B, Gauge=t.g, Check=false);
            z  = sv.z;
            nr = max(vecnorm(real(z).*t.g.e1 + imag(z).*t.g.e2, 2, 2));
            JA = reshape((  real(z).*t.g.e1 + imag(z).*t.g.e2 )', [], 1) / nr;
            JB = reshape(( -imag(z).*t.g.e1 + real(z).*t.g.e2 )', [], 1) / nr;
            switch kind
                case "rotating", D = (t.G*JA)*t.bank.a + (t.G*JB)*t.bank.b;
                case "standing", D = (t.G*JA)*t.bank.a;
            end
        end
    end

    methods (Test)

        function theBankHasTheShapeItClaims(t)
            n = t.bank.nAtoms;
            t.verifyEqual(n, 16*3);
            t.verifySize(t.bank.PA, [size(t.G,1) n]);
            t.verifySize(t.bank.PB, [size(t.G,1) n]);
            t.verifyEqual(numel(t.bank.vertex), n);
            t.verifyGreaterThan(min(t.bank.n0), 0);
            t.verifyEqual(t.bank.fc, 8*sqrt(2), 'RelTol', 0.02);
        end

        function aPlantedAtomIsFoundExactly(t)
            % the on-grid case: the atom IS in the dictionary, so this is a sanity check only
            k = 7;
            D = t.bank.PA(:,k)*t.bank.a + t.bank.PB(:,k)*t.bank.b;
            m = rheome.flow.vortexmatch(t.bank, D);
            t.verifyEqual(m.best.vertex, t.bank.vertex(k));
            t.verifyEqual(m.best.scaleMM, t.bank.scaleMM(k));
            t.verifyGreaterThan(m.best.energyFrac, 0.99);
        end

        function theMatchIsInvariantToTemporalPhase(t)
            % ⭐ the reason the energy needs no phase search
            k = 5;
            pA = t.bank.PA(:,k);  pB = t.bank.PB(:,k);
            e = zeros(1,4);  ph = [0 pi/3 pi/2 2.2];
            for i = 1:4
                aP = t.bank.a*cos(ph(i)) - t.bank.b*sin(ph(i));
                bP = t.bank.a*sin(ph(i)) + t.bank.b*cos(ph(i));
                m = rheome.flow.vortexmatch(t.bank, pA*aP + pB*bP);
                e(i) = m.best.energy;
                t.verifyEqual(m.best.vertex, t.bank.vertex(k));
            end
            t.verifyEqual(e, e(1)*ones(1,4), 'RelTol', 1e-6);
        end

        function noiseAloneScoresFarBelowAPlantedVortex(t)
            % ⭐ presence IS settled: measured ratio ~9299 at SNR 3
            D = t.i_plant(t.bank.vertex(4), 140, "rotating");
            mV = rheome.flow.vortexmatch(t.bank, D);
            rng(3);
            mN = rheome.flow.vortexmatch(t.bank, randn(size(D))*norm(D,'fro')/sqrt(numel(D)));
            t.verifyGreaterThan(mV.best.energyFrac, 0.5);
            t.verifyLessThan(mN.best.energyFrac, 0.05);
            t.verifyGreaterThan(mV.best.energyFrac/mN.best.energyFrac, 20);
        end

        function energyDoesNotSeparateRotatingFromStanding(t)
            % ⚠⚠ the caveat the docstring exists for. ⚠ It MUST be asserted over several
            %   locations: the per-location ratio runs 1.02 to 3.87, so a single seed proves
            %   nothing -- an earlier version of this test happened to pick the 3.87 one.
            locs = unique(t.bank.vertex, 'stable');
            locs = locs(1:min(8, numel(locs)));
            r = zeros(1, numel(locs));
            for i = 1:numel(locs)
                eR = rheome.flow.vortexmatch(t.bank, t.i_plant(locs(i), 140, "rotating")).best.energyFrac;
                eS = rheome.flow.vortexmatch(t.bank, t.i_plant(locs(i), 140, "standing")).best.energyFrac;
                r(i) = eR/eS;
            end
            t.verifyLessThan(median(r), 2, 'energy must not be read as evidence of rotation');
            t.verifyGreaterThan(max(r), 1, 'sanity: rotating should never score below standing');
        end

        function icohSeparatesThemAtEveryLocation(t)
            % ⭐ the discriminator that does reach rotation, asserted as a SEPARATION rather than
            %   a threshold: measured worst rotating 0.080 against best standing 0.023
            locs = unique(t.bank.vertex, 'stable');
            locs = locs(1:min(8, numel(locs)));
            iR = zeros(1, numel(locs));  iS = zeros(1, numel(locs));
            for i = 1:numel(locs)
                iR(i) = abs(rheome.flow.vortexmatch(t.bank, t.i_plant(locs(i), 140, "rotating")).best.icoh);
                iS(i) = abs(rheome.flow.vortexmatch(t.bank, t.i_plant(locs(i), 140, "standing")).best.icoh);
            end
            t.verifyGreaterThan(min(iR), max(iS), 'icoh must separate at EVERY location');
            t.verifyLessThan(max(iS), 0.1, 'a standing field must not read as rotating');
        end

        function mismatchedInputIsRejected(t)
            t.verifyError(@() rheome.flow.vortexmatch(t.bank, zeros(3, t.bank.nT)), ?MException);
            t.verifyError(@() rheome.flow.vortexmatch(t.bank, zeros(size(t.G,1), 7)), ?MException);
        end

    end
end

% Author: Diellor Basha, 2026
