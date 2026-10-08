classdef tFlowWindowTable < matlab.unittest.TestCase
% Invariants of rheome.flow.windowtable. The point of this table is that it admits NON-additive columns,
% so the tests pin down which properties survive that choice and which do not.
%
% Author: Diellor Basha, 2026
    properties
        T; X
    end
    methods (TestClassSetup)
        function build(t)
            try
                B = rheome.load.bases(rheomeTestSubject()); %#ok<NASGU>
            catch
                t.assumeFail('cached Brainstorm data for test subject is not present');
            end
            % ⚠ MaxDepth must reach 4: depths 1-3 are 133 mm and up, so a tree stopped at 3 has no
            % sub-resolution tiles at all and the admissibility test has nothing to find.
            [t.T, t.X] = rheome.flow.windowtable(rheomeTestSubject(), WindowSec=2, MaxDepth=4, ...
                Bands=[8 16], MaxWindows=3, Flow=false, Verbose=false);
        end
    end
    methods (Test)
        function energyRollsUpExactly(t)
            % ⭐ the ONE column for which a parent is the sum of its children
            for w = unique(t.T.window)'
                m = t.T.window == w;
                e1 = sum(t.T.energy(m & t.T.depth==1));
                e3 = sum(t.T.energy(m & t.T.depth==max(t.T.depth)));
                t.verifyEqual(e1, e3, 'RelTol', 1e-9, ...
                    'rheome.geom.tree is a disjoint partition, so energy must roll up');
            end
        end
        function perTileColumnsVaryByTile(t)
            % ⚠ this is the regression test for a real bug: the gauge split and the vortex count
            % were computed hemisphere-wide and copied into every row, where they looked per-tile.
            m = t.T.window==1 & t.T.depth==3;
            for c = ["meanAct" "peakAct" "crest" "energy" "normalShare" "tangShare"]
                t.verifyGreaterThan(numel(unique(t.T.(c)(m))), 1, ...
                    sprintf('%s is replicated across tiles, so it is not a per-tile measurement', c));
            end
        end
        function windowPrefixedColumnsAreConstantWithinAWindow(t)
            m = t.T.window==1;
            for c = string(t.T.Properties.VariableNames)
                if startsWith(c, "win_")
                    t.verifyEqual(numel(unique(t.T.(c)(m))), 1, ...
                        sprintf('%s carries the win_ prefix but varies by tile', c));
                end
            end
        end
        function theGaugeSharesArePartitions(t)
            s = t.T.normalShare + t.T.tangShare;
            t.verifyEqual(s, ones(size(s)), 'AbsTol', 1e-9);
        end
        function admissibilityIsTheResolutionFloorAndRowsAreKept(t)
            % ⚠ inadmissible tiles must be PRESENT and flagged, never dropped: they are the control
            t.verifyEqual(t.T.admissible, t.T.nodeDiameterMM > 104);
            t.verifyTrue(any(~t.T.admissible), 'the sub-resolution control rows were dropped');
        end
        function eigenmodeGroupsAreBlocksOfFour(t)
            t.verifyEqual(numel(t.X.wavelengthMM), size(t.X.GroupPower,1));
            t.verifyEqual(mod(numel(t.X.octave), 1), 0);
            t.verifyTrue(all(t.X.octave >= 1 & t.X.octave <= 6));
        end
        function everyRowIsOneWindowBandTile(t)
            k = unique(t.T(:, {'window','bandIdx','node'}));
            t.verifyEqual(height(k), height(t.T), 'the index is not unique');
        end
    end
end
