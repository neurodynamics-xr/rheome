classdef tSensorTree < matlab.unittest.TestCase
% The position pyramid: a binary tree whose leaves are the sensors, whose internal nodes
% partition their members between two children, and whose cuts on a lattice follow the
% lattice's axes.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function leavesAreTheSensorsAndNodesPartitionTheirChildren(tc)
            arr = rheome.sensors.grid('Size', [6 8], 'Pitch', 1e-2);
            G = rheome.sensors.graph(arr, 'Faces', false);
            T = rheome.sensors.tree(G);
            tc.verifyEqual(nnz(T.is_leaf), arr.nCh);
            tc.verifyEqual(sort(T.channel_id(T.is_leaf))', 1:arr.nCh);
            tc.verifyEqual(T.parent_id(1), 0);  tc.verifyEqual(T.depth(1), 0);
            tc.verifyEqual(T.members{1}, 1:arr.nCh);
            for i = find(~T.is_leaf)'
                kids = find(T.parent_id == T.node_id(i));
                tc.verifyNumElements(kids, 2, sprintf('node %d', i));
                tc.verifyEqual(sort([T.members{kids(1)}, T.members{kids(2)}]), sort(T.members{i}));
                tc.verifyTrue(isempty(intersect(T.members{kids(1)}, T.members{kids(2)})));
                tc.verifyTrue(all(T.depth(kids) == T.depth(i) + 1));
                tc.verifyTrue(all(kids > i));                                    % children after parents
            end
            tc.verifyEqual(height(T), 2 * arr.nCh - 1);                          % a full binary tree
        end

        function theFirstCutOfARectangleFollowsItsLongAxis(tc)
            arr = rheome.sensors.grid('Size', [4 12], 'Pitch', 1e-2);                   % long along y
            G = rheome.sensors.graph(arr, 'Faces', false);
            T = rheome.sensors.tree(G);
            tc.verifyEqual(T.axis{1}, 'y');
            kids = find(T.parent_id == 1);
            tc.verifyEqual(sort(T.n_sensors(kids))', [24 24]);
        end

        function minSizeStopsTheSplit(tc)
            arr = rheome.sensors.grid('Size', [4 4], 'Pitch', 1e-2);
            T = rheome.sensors.tree(rheome.sensors.graph(arr, 'Faces', false), MinSize=4);
            tc.verifyTrue(all(T.n_sensors(T.is_leaf) <= 4));
            tc.verifyTrue(all(T.n_sensors(~T.is_leaf) > 4));
        end

    end
end
