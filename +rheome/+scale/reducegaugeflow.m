function [Gr, Pl] = reducegaugeflow(inDir, outDir)
% SCALE.REDUCEGAUGEFLOW  Per-subject gaugeflow.csv -> the group alpha flow in the shared frame.
%
%   [Gr, Pl] = rheome.scale.reducegaugeflow(inDir, outDir)
%
% Reads every *gaugeflow.csv under inDir (rheome.scale.run's per-subject folders, or nsp's flat
% <sub>_gaugeflow.csv) and, per cohort x band x ico3 patch, pools SUBJECTS -- each subject's patch
% mean counts once. Each patch is read in chart "z" (poles at FreeSurfer's sphere +-z, the group
% gauge) unless z excludes it (rheome.geom.spherepatches: a pole, or the frame turning > 15 deg across it),
% in which case chart "x" is used and .chart says so. DESCRIPTIVE: no test, no correction.
%
% Writes outDir/group_gaugeflow.csv (Gr), one row per cohort x band x patch:
%   n_subjects
%   vn, vw            group mean of the subjects' mean north / west velocity, m/s
%   heading_deg       atan2(vw, vn): 0 = north (towards the +pole), 90 = west
%   speed_of_mean     |(vn, vw)|, m/s       mean_speed   mean over subjects of |(vn_s, vw_s)|, m/s
%   R                 speed_of_mean / mean_speed: 1 = every subject's mean flow points one way
%   axis_deg          principal axis of the group mean orientation tensor, in [-90 90), 0 = north
%   anisotropy        (l1 - l2)/(l1 + l2) of that tensor: 0 = no preferred axis
%   env_median        median over subjects of the patch's mean activation
%   colat_deg lon_deg of the patch in the chart used
% and outDir/poles.csv (Pl), one row per subject x z-pole patch: the patch's env_mean percentile
% among that subject's patches (periodic band and total) -- how active the cortex under each
% group-gauge pole is.
%
% ⚠ CONP AND PREVENT-AD ARE THE SAME PEOPLE (rheome.scale.reduce): cohorts are reduced separately and
% never pooled here.
%
% See also: rheome.scale.run, rheome.scale.measure_periodicflow, rheome.geom.sphereframe, rheome.geom.spherepatches, rheome.scale.reduce
%
% Author: Diellor Basha, 2026

    arguments
        inDir (1,:) char
        outDir (1,:) char
    end
    f = dir(fullfile(inDir, '**', '*gaugeflow.csv'));
    if isempty(f), error('scale:reducegaugeflow:none', 'no *gaugeflow.csv under %s', inDir); end
    D = table();
    for k = 1:numel(f)
        t = readtable(fullfile(f(k).folder, f(k).name), 'TextType', 'string');
        D = [D; t]; %#ok<AGROW>
    end
    D.cohort(ismissing(D.cohort) | D.cohort == "") = "unlabelled";
    % the shared frame: z where it holds, x where z is excluded (same patch numbering in both)
    Z = D(D.chart == "z", :);  Xc = D(D.chart == "x", :);
    useX = Z.excluded == 1;
    [~, ix] = ismember(Z(useX, {'subject','band','patch'}), Xc(:, {'subject','band','patch'}));
    S = [Z(~useX, :); Xc(ix(ix > 0), :)];
    S = S(S.n_samples > 0 & isfinite(S.vn_mean), :);

    [g, cohort, band, patch] = findgroups(S.cohort, S.band, S.patch);
    sp = sqrt(S.vn_mean.^2 + S.vw_mean.^2);
    n  = splitapply(@numel, S.vn_mean, g);
    vn = splitapply(@mean, S.vn_mean, g);  vw = splitapply(@mean, S.vw_mean, g);
    ms = splitapply(@mean, sp, g);
    tnn = splitapply(@mean, S.tnn, g);  tnw = splitapply(@mean, S.tnw, g);  tww = splitapply(@mean, S.tww, g);
    env = splitapply(@median, S.env_mean, g);
    chart = splitapply(@(c) c(1), S.chart, g);
    colat = splitapply(@(c) c(1), S.colat_deg, g);  lon = splitapply(@(c) c(1), S.lon_deg, g);
    l = sqrt(((tnn - tww)/2).^2 + tnw.^2);  tr = (tnn + tww)/2;
    som = hypot(vn, vw);
    Gr = table(cohort, band, patch, chart, n, vn, vw, atan2d(vw, vn), som, ms, som ./ ms, ...
               0.5*atan2d(2*tnw, tnn - tww), l ./ tr, env, colat, lon, 'VariableNames', ...
               {'cohort','band','patch','chart','n_subjects','vn','vw','heading_deg','speed_of_mean', ...
                'mean_speed','R','axis_deg','anisotropy','env_median','colat_deg','lon_deg'});

    % how active is the cortex under each z pole, per subject
    Pl = table();
    for s = unique(Z.subject)'
        for b = unique(Z.band)'
            z = Z(Z.subject == s & Z.band == b, :);
            for r = find(z.has_pole == 1)'
                Pl = [Pl; table(s, z.cohort(r), b, z.patch(r), z.colat_deg(r), ...
                      100*mean(z.env_mean < z.env_mean(r), 'omitnan'), 'VariableNames', ...
                      {'subject','cohort','band','patch','colat_deg','env_percentile'})]; %#ok<AGROW>
            end
        end
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    writetable(Gr, fullfile(outDir, 'group_gaugeflow.csv'));
    writetable(Pl, fullfile(outDir, 'poles.csv'));
end

% Author: Diellor Basha, 2026
