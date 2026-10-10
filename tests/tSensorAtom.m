classdef tSensorAtom < matlab.unittest.TestCase
    % rheome.forward.sensoratom: each sensor's sensitivity fitted as a graph wavelet.
    %
    % Author: Diellor Basha, 2026

    properties
        G; db; S
    end

    methods (TestClassSetup)
        function build(t)
            rheomeTestSubject(t, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try
                B  = rheome.load.bases(rheomeTestSubject());
                st = rheome.load.study(rheomeTestSubject());
                d  = rheome.load.dirac(rheomeTestSubject());
            catch
                t.assumeFail('cached bases/study/dirac for test subject are not present');
            end
            gv   = double(B.L.gv(:))';
            t.S  = B.L.S;
            t.G  = rheome.forward.leadfield(st, GlobalVertices=gv);
            cols = find(d.Hemisphere == 1);
            t.db = struct('Phi', d.Phi(reshape((gv'-1)*4+(1:4),[],1), cols), ...
                          'Lambda', d.Lambda(cols), 'nVert', size(t.S.Vertices,1), ...
                          'nModes', numel(cols));
        end
    end

    methods (Test)

        function itFitsEverySensor(t)
            o = rheome.forward.sensoratom(t.G, t.db, struct('Surface',t.S,'Candidates',40));
            nCh = size(t.G,1);
            t.verifyEqual(height(o.table), nCh);
            t.verifySize(o.varExplained, [nCh 1]);
            t.verifyTrue(all(o.varExplained >= 0 & o.varExplained <= 1+1e-9));
            t.verifyTrue(all(ismember(o.vertex(:,1), o.candidates)));
        end

        function aWaveletBeatsAVortexByFarOnTheSameRows(t)
            % ⭐ the measurement the function exists for: 0.69 against a vortex's 0.087
            o = rheome.forward.sensoratom(t.G, t.db, struct('Surface',t.S,'Candidates',60));
            t.verifyGreaterThan(median(o.varExplained), 0.4, ...
                'a graph wavelet should explain far more of a row than a vortex does');
        end

        function moreAtomsExplainMore(t)
            o1 = rheome.forward.sensoratom(t.G, t.db, struct('Surface',t.S,'Candidates',40,'NAtoms',1));
            o3 = rheome.forward.sensoratom(t.G, t.db, struct('Surface',t.S,'Candidates',40,'NAtoms',3));
            t.verifySize(o3.vertex, [size(t.G,1) 3]);
            t.verifyGreaterThanOrEqual(median(o3.varExplained), median(o1.varExplained));
        end

        function explicitCandidatesAreHonoured(t)
            cl = [100 500 1000 2000 4000];
            o  = rheome.forward.sensoratom(t.G, t.db, struct('Candidates', cl));
            t.verifyEqual(sort(o.candidates), sort(cl));
            t.verifyTrue(all(ismember(o.vertex(:,1), cl)));
        end

        function badInputsAreRejected(t)
            t.verifyError(@() rheome.forward.sensoratom(t.G, t.db, ...
                struct('Family','gaussian')), ?MException);
            t.verifyError(@() rheome.forward.sensoratom(t.G, t.db, ...
                struct('Candidates', 20)), ?MException);   % a count with no Surface
        end

    end
end

% Author: Diellor Basha, 2026
