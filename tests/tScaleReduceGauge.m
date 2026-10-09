classdef tScaleReduceGauge < matlab.unittest.TestCase
% rheome.scale.reducegauge on hand-made gaugetensor.csv files: the group tensor's axis, the axis agreement
% across subjects, the polar exclusion of the axis, and one row per dataset.
%
% Author: Diellor Basha, 2026

    methods (Test)
        function groupTensorAndPolarExclusion(tc)
            d = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            % three subjects in dataset omega share tile 1's axis at 30 deg; tile 2 is polar
            ax = [30 30 30];  cl = [80 25];  coh = ["omega" "omega" "omega" "preventad-meg"];
            for s = 1:4
                p = fullfile(d, 'cf-atlas', sprintf('sub-%02d', s));  mkdir(p);  pre = sprintf('sub-%02d_', s);
                if s == 1, pre = ''; end                    % one in run's layout, the rest in nsp's flat one
                a = [ax 120];  a = a(s);  an = 0.5;          % tensor = (I + an*[cos2a sin2a; sin2a -cos2a])/2
                t11 = (1 + an*cosd(2*a))/2;  t22 = (1 - an*cosd(2*a))/2;  t12 = an*sind(2*a)/2;
                t = table(["L";"L"], [3;3], [1;2], [10;10], [t11;t11], [t22;t22], [t12;t12], [an;an], [a;a], ...
                          [0.2;0.2], [1;1], cl(:), 'VariableNames', {'hemi','depth','node_id','n_vertices', ...
                          't11','t22','t12','anisotropy','axis_deg','normal_share','trace','colat_median_deg'});
                writetable(t, fullfile(p, [pre 'gaugetensor.csv']));
                writetable(table(string(p), coh(s), 'VariableNames', {'subject','dataset'}), fullfile(p, [pre 'metrics.csv']));
            end
            G = rheome.scale.reducegauge(d, fullfile(d, 'out'));
            tc.verifyTrue(isfile(fullfile(d, 'out', 'gaugegroup.csv')));
            r = G(G.dataset == "omega" & G.node_id == 1, :);
            tc.verifyEqual(r.n, 3);
            tc.verifyEqual(r.axis_deg, 30, 'AbsTol', 1e-9);
            tc.verifyEqual(r.axis_R, 1, 'AbsTol', 1e-12);
            tc.verifyEqual(r.anisotropy, 0.5, 'AbsTol', 1e-12);
            tc.verifyFalse(r.polar);
            q = G(G.dataset == "omega" & G.node_id == 2, :);
            tc.verifyTrue(q.polar);  tc.verifyTrue(isnan(q.axis_deg));  tc.verifyEqual(q.anisotropy, 0.5, 'AbsTol', 1e-12);
            m = G(G.dataset == "preventad-meg" & G.node_id == 1, :);
            tc.verifyEqual(m.axis_deg, -60, 'AbsTol', 1e-9);  % 120 deg is the same axis as -60
        end
    end
end

% Author: Diellor Basha, 2026
