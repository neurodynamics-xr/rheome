classdef tSensorBank < matlab.unittest.TestCase
    % rheome.forward.sensorbank: a Helmholtz graph-wavelet bank fitted to leadfield rows.
    %
    % Author: Diellor Basha, 2026

    properties
        G; S; lbo; sens
    end

    methods (TestClassSetup)
        function build(t)
            try
                B  = rheome.load.bases(rheomeTestSubject());
                st = rheome.load.study(rheomeTestSubject());
            catch
                t.assumeFail('cached bases/study for test subject are not present');
            end
            t.S   = B.L.S;
            t.lbo = B.L.lbo;
            t.G   = rheome.forward.leadfield(st, GlobalVertices=double(B.L.gv(:))');
            rng(7);  t.sens = randperm(size(t.G,1), 4);
        end
    end

    methods (Test)

        function theOutputHasTheShapeItClaims(t)
            o = rheome.forward.sensorbank(t.G, t.S, t.lbo, ...
                    struct('NAtoms',3,'Candidates',60,'Sensors',t.sens));
            t.verifyEqual(height(o.atoms), 3*numel(t.sens));
            t.verifyTrue(all(ismember(o.atoms.family, ["grad","rot","norm"])));
            t.verifyTrue(all(ismember(o.atoms.vertex, o.candidates)));
            t.verifyTrue(all(o.ve(t.sens) > 0 & o.ve(t.sens) <= 1+1e-9));
        end

        function moreAtomsExplainMore(t)
            % ⭐ orthogonal matching pursuit must be monotone in the field domain
            o1 = rheome.forward.sensorbank(t.G, t.S, t.lbo, ...
                    struct('NAtoms',1,'Candidates',60,'Sensors',t.sens));
            o5 = rheome.forward.sensorbank(t.G, t.S, t.lbo, ...
                    struct('NAtoms',5,'Candidates',60,'Sensors',t.sens));
            t.verifyGreaterThan(median(o5.ve(t.sens)), median(o1.ve(t.sens)));
        end

        function allThreeFamiliesGetUsed(t)
            % ⚠ a leadfield row is 0.41 gradient, 0.16 curl, 0.34 normal, so a fit that
            %   never selects the rotor family has lost the circulation by construction
            o = rheome.forward.sensorbank(t.G, t.S, t.lbo, ...
                    struct('NAtoms',10,'Candidates',60,'Sensors',t.sens));
            t.verifyGreaterThan(sum(o.atoms.family=="rot"), 0, 'no rotor atom was ever chosen');
            t.verifyGreaterThan(numel(unique(o.atoms.family)), 1);
        end

        function itBeatsTheVectorAtomItReplaces(t)
            % §41: the Dirac vector atom scored ve 0.104 and curl r 0.110 on these same rows
            o = rheome.forward.sensorbank(t.G, t.S, t.lbo, ...
                    struct('NAtoms',10,'Candidates',60,'Sensors',t.sens));
            t.verifyGreaterThan(median(o.ve(t.sens)), 0.104);
            t.verifyGreaterThan(median(o.curlCorr(t.sens)), 0.110);
        end

        function theReconstructionMatchesTheReportedFit(t)
            o = rheome.forward.sensorbank(t.G, t.S, t.lbo, ...
                    struct('NAtoms',5,'Candidates',60,'Sensors',t.sens,'Reconstruct',true));
            t.verifySize(o.Bhat, size(t.G));
            for c = t.sens
                L = t.G(c,:)';
                veDirect = 1 - norm(L - o.Bhat(c,:)')^2/norm(L)^2;
                t.verifyEqual(veDirect, o.ve(c), 'AbsTol', 1e-8);
            end
        end

    end
end

% Author: Diellor Basha, 2026
