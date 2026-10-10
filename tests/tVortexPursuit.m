classdef tVortexPursuit < matlab.unittest.TestCase
    % rheome.flow.vortexpursuit: matching pursuit over the vortex dictionary.
    %
    % Author: Diellor Basha, 2026

    properties
        B; g; st; bank; G; fs
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
            t.fs = t.bank.fs;
        end
    end

    methods (Access = private)
        function [D, v0] = i_plant(t, nT)
            v0 = t.bank.vertex(4);
            k  = find(t.bank.vertex == v0 & t.bank.scaleMM == 140, 1);
            tt = (-fix(nT/2):fix((nT-1)/2))/t.fs;
            fc = t.bank.fc;
            psi = ifftshift(exp(2i*pi*fc*tt).*exp(-tt.^2/(2*(15/(2*pi*fc))^2)));
            psi = psi/norm(psi);
            D = zeros(size(t.G,1), nT);  rng(1);
            for c0 = 1:round(1.5*t.fs):nT
                ps = circshift(psi, c0-1)*exp(2i*pi*rand);
                D  = D + t.bank.PA(:,k)*real(ps) + t.bank.PB(:,k)*imag(ps);
            end
        end
    end

    methods (Test)

        function theOutputHasTheShapeItClaims(t)
            nT = round(6*t.fs);
            pu = rheome.flow.vortexpursuit(t.bank, t.i_plant(nT), 'NAtoms', 5);
            t.verifyEqual(height(pu.atoms), 5);
            t.verifySize(pu.Bhat, [size(t.G,1) nT]);
            t.verifyEqual(numel(pu.resid), 5);
            t.verifyTrue(all(ismember(pu.atoms.cycles, [2 5 15])));
        end

        function theResidualFallsMonotonically(t)
            % ⭐ the defining property of matching pursuit: every atom reduces total error
            nT = round(6*t.fs);
            pu = rheome.flow.vortexpursuit(t.bank, t.i_plant(nT), 'NAtoms', 8);
            t.verifyTrue(all(diff(pu.resid) <= 1e-12), 'residual must not increase');
            t.verifyLessThan(pu.resid(end), pu.resid(1));
            t.verifyLessThanOrEqual(pu.resid(1), 1);
        end

        function aPlantedAtomIsFoundFirst(t)
            nT = round(6*t.fs);
            [D, v0] = t.i_plant(nT);
            pu = rheome.flow.vortexpursuit(t.bank, D, 'NAtoms', 3);
            t.verifyEqual(pu.atoms.vertex(1), v0);
            t.verifyEqual(pu.atoms.scaleMM(1), 140);
        end

        function heldOutChannelsDoNotEnterTheFit(t)
            % ⚠ the whole point of the benchmark: the fit must be blind to those rows
            nT = round(6*t.fs);
            D  = t.i_plant(nT);
            ch = 1:2:size(D,1);
            D2 = D;  D2(setdiff(1:size(D,1), ch), :) = 0;   % wreck the held-out rows
            p1 = rheome.flow.vortexpursuit(t.bank, D,  'NAtoms', 4, 'Channels', ch);
            p2 = rheome.flow.vortexpursuit(t.bank, D2, 'NAtoms', 4, 'Channels', ch);
            t.verifyEqual(p1.atoms.vertex, p2.atoms.vertex, 'the fit used held-out channels');
            t.verifyEqual(p1.atoms.c0, p2.atoms.c0, 'RelTol', 1e-10);
        end

        function theSourceEstimateForwardsBackToTheSensorFit(t)
            % ⭐ J is what makes it comparable with a distributed inverse at all
            nT = round(4*t.fs);
            pu = rheome.flow.vortexpursuit(t.bank, t.i_plant(nT), 'NAtoms', 3, ...
                     'Source', true, 'Name', string(rheomeTestSubject()), 'Bases', t.B, 'Gauge', t.g);
            t.verifySize(pu.J, [3*size(t.B.L.S.Vertices,1) nT]);
            t.verifyEqual(t.G*pu.J, pu.Bhat, 'RelTol', 1e-8, 'AbsTol', 1e-20);
        end

        function sourceNeedsAName(t)
            nT = round(4*t.fs);
            t.verifyError(@() rheome.flow.vortexpursuit(t.bank, t.i_plant(nT), ...
                'NAtoms', 1, 'Source', true), ?MException);
        end

    end
end

% Author: Diellor Basha, 2026
