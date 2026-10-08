classdef tScaleReduce < matlab.unittest.TestCase
% TSCALEREDUCE  rheome.scale.rows and rheome.scale.reduce on synthetic per-subject folders (no cached data).
%
% Author: Diellor Basha, 2026

    methods (Test)
        function rowsShape(tc)
            T = rheome.scale.rows("resolution", ["a" "b"], [1 2], "mm");
            tc.verifyEqual(height(T), 2);
            tc.verifyEqual(T.band, ["";""]);
            tc.verifyEqual(T.value, [1;2]);
        end

        function reduceDistributionsAnchorAndTwins(tc)
            d = tempname;  tc.addTeardown(@() rmdir(d, 's'));
            nm = ["a_sub1" "a_sub2" "a_sub3" "a_subB1" "a_subB2" "a_subB3" ...
                  "b_subT1" "b_subT2" "b_subT3" "a_ref"];
            co = ["norm" "norm" "norm" "B" "B" "B" "T" "T" "T" "norm"];
            tw = ["" "" "" "sub-T1" "sub-T2" "sub-T3" "sub-B1" "sub-B2" "sub-B3" ""];
            v  = [10 20 30 1 2 3 1 2 3 25];
            for i = 1:numel(nm)
                p = fullfile(d, nm(i));  mkdir(p);
                M = [table(nm(i), co(i), "x", 'VariableNames', {'subject','cohort','dataset'}) rheome.scale.rows("res", "r50", v(i), "mm")];
                writetable(M, fullfile(p, 'metrics.csv'));
                writetable(table("res", 1, 1, "ok", "", 'VariableNames', {'analysis','seconds','peak_rss_GB','status','message'}), fullfile(p, 'timing.csv'));
                writelines("twin" + char(9) + tw(i), fullfile(p, 'job.tsv'));
            end
            G = rheome.scale.reduce(d, fullfile(d, '_out'), Anchor="a_ref", TwinSide="a_");
            r = G.distribution(G.distribution.cohort == "norm", :);
            tc.verifyEqual(r.n, 3);                        % the anchor is not counted in its reference
            tc.verifyEqual(r.median, 20);
            tc.verifyEqual(G.anchor.percentile, 200/3, 'AbsTol', 1e-9);
            tc.verifyEqual(G.twins.n_pairs, 3);            % each pair once, from the TwinSide only
            tc.verifyEqual(G.twins.spearman, 1, 'AbsTol', 1e-12);
            tc.verifyEqual(G.status.ok, 10);
            tc.verifyFalse(any(G.cohorts.cohort == "anchor"));
        end
    end
end

% Author: Diellor Basha, 2026
