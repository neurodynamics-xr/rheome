classdef tVortexTrack < matlab.unittest.TestCase
    % rheome.flow.vortextrack: vortex parameters as functions of time.
    %
    % Author: Diellor Basha, 2026

    properties
        B; g; st; bank; G; fs; fc
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
            t.G    = rheome.forward.leadfield(t.st, GlobalVertices=double(t.B.L.gv(:))');
            t.bank = rheome.flow.vortexbank(rheomeTestSubject(), Locations=20, Bases=t.B, Gauge=t.g, ...
                                     Study=t.st, Verbose=false);
            t.fs = t.bank.fs;  t.fc = t.bank.fc;
        end
    end

    methods (Access = private)
        function D = i_plant(t, vtx, lam, ncyc, nT)
            sv = rheome.flow.seedvortex(rheomeTestSubject(), Vertex=vtx, WavelengthMM=lam, ...
                                 Bases=t.B, Gauge=t.g, Check=false);
            z  = sv.z;
            nr = max(vecnorm(real(z).*t.g.e1 + imag(z).*t.g.e2, 2, 2));
            JA = reshape((  real(z).*t.g.e1 + imag(z).*t.g.e2 )', [], 1)/nr;
            JB = reshape(( -imag(z).*t.g.e1 + real(z).*t.g.e2 )', [], 1)/nr;
            tt  = (-fix(nT/2):fix((nT-1)/2))/t.fs;
            sd  = ncyc/(2*pi*t.fc);
            psi = ifftshift(exp(2i*pi*t.fc*tt).*exp(-tt.^2/(2*sd^2)));  psi = psi/norm(psi);
            e = abs(psi).^2; [~,o] = sort(e,'descend');
            sp = find(cumsum(e(o)) >= 0.999*sum(e), 1);
            D = zeros(size(t.G,1), nT);  rng(2);
            for c0 = 1:max(1,round(sp/2)):nT
                ps = circshift(psi, c0-1)*exp(2i*pi*rand);
                D  = D + (t.G*JA)*real(ps) + (t.G*JB)*imag(ps);
            end
        end
    end

    methods (Test)

        function theTrackHasTheShapeItClaims(t)
            nT = round(8*t.fs);
            D  = t.i_plant(t.bank.vertex(3), 140, 15, nT);
            tr = rheome.flow.vortextrack(t.bank, D, 'Stride', 0.5);
            nE = numel(tr.t);
            t.verifyEqual(height(tr.table), nE);
            t.verifyEqual(numel(tr.vertex), nE);
            t.verifyTrue(all(ismember(tr.cycles, [2 5 15])));
            t.verifyTrue(all(tr.energyFrac >= 0));
            t.verifyEqual(tr.t(2)-tr.t(1), 0.5, 'RelTol', 0.05);
        end

        function aPlantedSourceIsLocated(t)
            nT = round(8*t.fs);
            v0 = t.bank.vertex(3);
            tr = rheome.flow.vortextrack(t.bank, t.i_plant(v0, 140, 15, nT), 'Stride', 0.5, 'Icoh', false);
            k  = tr.t > 1 & tr.t < 7;
            t.verifyEqual(mode(tr.vertex(k)), v0);
            t.verifyEqual(mode(tr.scaleMM(k)), 140);
        end

        function silenceScoresFarBelowAnActiveStretch(t)
            nT = round(8*t.fs);
            D  = t.i_plant(t.bank.vertex(3), 140, 15, nT);
            trA = rheome.flow.vortextrack(t.bank, D, 'Stride', 0.5, 'Icoh', false);
            trS = rheome.flow.vortextrack(t.bank, zeros(size(D))+1e-18*randn(size(D)), ...
                                   'Stride', 0.5, 'Icoh', false);
            t.verifyGreaterThan(median(trA.energyFrac), 10*median(trS.energyFrac));
        end

        function theLengthPenaltyIsWhatMakesQIdentifiable(t)
            % ⚠⚠ measured: raw energy returns the LONGEST atom whatever the truth, 0/6 correct for
            %    planted Q of 2 and 5 at SNR 3. This pins that, so the default cannot drift back.
            nT = round(10*t.fs);
            D  = t.i_plant(t.bank.vertex(3), 140, 2, nT);
            k  = @(tr) tr.t > 1.5 & tr.t < 8.5;
            trP = rheome.flow.vortextrack(t.bank, D, 'Stride', 0.25, 'Icoh', false);              % 0.5
            trR = rheome.flow.vortextrack(t.bank, D, 'Stride', 0.25, 'Icoh', false, 'Penalty', 0);
            t.verifyEqual(mode(trP.cycles(k(trP))), 2, 'the penalty should recover a short atom');
            t.verifyEqual(mode(trR.cycles(k(trR))), 15, 'raw energy should pin to the longest');
        end

        function aLongPlantedAtomIsRecoveredToo(t)
            % the penalty must not simply always answer "short"
            nT = round(10*t.fs);
            tr = rheome.flow.vortextrack(t.bank, t.i_plant(t.bank.vertex(3), 140, 15, nT), ...
                                  'Stride', 0.25, 'Icoh', false);
            k = tr.t > 1.5 & tr.t < 8.5;
            t.verifyEqual(mode(tr.cycles(k)), 15);
        end

        function mismatchedInputIsRejected(t)
            t.verifyError(@() rheome.flow.vortextrack(t.bank, zeros(3, 100)), ?MException);
        end

    end
end

% Author: Diellor Basha, 2026
