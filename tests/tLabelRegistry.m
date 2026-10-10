classdef tLabelRegistry < matlab.unittest.TestCase
% The calibration registry (rheome.select.labelkinds), its guard in rheome.select.measure, and the whole-cortex
% node ids (rheome.geom.cortexnodes) every cortical label is addressed by.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function everyKindHasAKnownCalibrationAndEvidence(tc)
            K = rheome.select.labelkinds();
            tc.verifyTrue(all(ismember(K.calibration, ["validated" "compressed" "pending" "unmeasurable"])));
            tc.verifyTrue(all(strlength(K.evidence) > 0));
            tc.verifyEqual(numel(unique(K.kind)), height(K));
            tc.verifyEqual(rheome.select.labelkinds("blob_size").calibration, "unmeasurable");
            tc.verifyEmpty(rheome.select.labelkinds("no_such_kind"));
        end

        function measureRefusesAnUnmeasurableKind(tc)
            db = struct('recording_id', "x", 'grid', struct('Lmax', 2, 'K', [4 2 1]), 'meta', struct('C', 1), ...
                        'groupNodes', [], 'labelFile', [tempname '.mat']);
            rows = table(0, 1, 1, 150, 'VariableNames', {'level','k','unit_id','value'});
            tc.verifyError(@() rheome.select.measure(db, rows, Kind="envelope_scale", Scope="cortex"), ...
                'select:measure:unmeasurable');
        end

        function cortexIdsAreAWholeCortexHeap(tc)
            rheomeTestSubject(tc, 'subject');   % skips, with the reason, when the cache or the subject is absent
            try, C = rheome.geom.cortexnodes(string(rheomeTestSubject()), MaxDepth=4);
            catch, tc.assumeFail('cached bases for the test subject are not present'); end
            N = C.nodes;
            tc.verifyEqual(numel(unique(N.node_id)), height(N));
            tc.verifyEqual(N.node_id(N.depth == 1)', [2 3]);                 % the two hemispheres
            nr = N.node_id > 1;
            tc.verifyEqual(N.parent_id(nr), floor(N.node_id(nr) / 2));
            tc.verifyEqual(N.depth(nr), floor(log2(N.node_id(nr))));        % depth by arithmetic
            tc.verifyTrue(all(N.hemi(N.node_id >= 2 & N.node_id < 3*2.^floor(log2(N.node_id)) ... 
                & floor(N.node_id ./ 2.^(floor(log2(N.node_id))-1)) == 2) == "L"));
            tc.verifyEqual(C.id(3, 1), 3);  tc.verifyEqual(C.id(2, 5), 9);  tc.verifyEqual(C.id(3, 5), 13);
            % area adds up: each internal node is the sum of its two children
            for i = find(~N.is_leaf)'
                kids = N.parent_id == N.node_id(i);
                tc.verifyEqual(sum(N.area(kids)), N.area(i), 'RelTol', 1e-10);
            end
        end
    end
end

% Author: Diellor Basha, 2026
