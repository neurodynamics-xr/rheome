function G = reducegauge(inDir, outDir, opts)
% SCALE.REDUCEGAUGE  Per-subject gaugetensor.csv -> the group alpha orientation tensor per shared tile.
%
%   G = rheome.scale.reducegauge(inDir, outDir)
%   G = rheome.scale.reducegauge(inDir, outDir, PolarDeg=35)
%
% inDir holds rheome.scale.run's per-subject folders (<sub>/gaugetensor.csv, <sub>/metrics.csv) or nsp cf-atlas's
% flat <sub>_gaugetensor.csv and <sub>_metrics.csv (any depth); metrics.csv gives the dataset (OMEGA, PREVENT-AD),
% which is the group: cf-atlas labels every subject's cohort "unlabelled". The tiles are the template's
% group tiles reached through the registered sphere, and the gauge is rheome.geom.sphereframe with the poles at
% the FreeSurfer sphere's +-z in every participant (MS1 section 9, 21cc8d9b), so tile k and its north are the
% same in everyone and the subjects' tensors can be averaged directly. Writes outDir/gaugegroup.csv, one row per
% (dataset, hemi, depth, node_id):
%   n               subjects with the tile
%   t11 t22 t12     the mean of the subjects' trace-normalised tensors (each subject weighs the same)
%   anisotropy      (l1-l2)/(l1+l2) of that mean tensor; axis_deg its major axis (0 = north, 90 = west)
%   axis_R          axis agreement across subjects: |mean exp(2i axis)| (1 = all subjects share the axis)
%   anisotropy_subject_median, normal_share_median, colat_median_deg   medians over subjects
%   polar           the tile's median colatitude is within PolarDeg of a pole
%
% ⚠ A POLE IN ACTIVE CORTEX IS EXCLUDED FROM THE AXIS, NOT FROM THE TILE. The +z pole is precentral, in active
% cortex; near it the meridian turns through 360 deg inside one tile, so the axis there is the frame's, not the
% current's. For polar tiles axis_deg and axis_R are NaN. Anisotropy and normal share are gauge-invariant and are
% kept. ⚠ THE THRESHOLD IS ON THE MEDIAN. A depth-3 tile is 1/8 of the sphere (pi/2 sr, a 41 deg cap), so a tile
% that holds a pole has median colatitude <= ~29 deg; PolarDeg 35 flags it, a tile beside the cap is not flagged.
% Descriptive only: no test across tiles or cohorts.
%
% See also: rheome.scale.measure_atlas, rheome.geom.sphereframe, rheome.scale.reduce
%
% Author: Diellor Basha, 2026

    arguments
        inDir (1,:) char
        outDir (1,:) char
        opts.PolarDeg (1,1) double = 35
    end
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    f = dir(fullfile(inDir, '**', '*gaugetensor.csv'));      % <sub>/gaugetensor.csv (run) or <sub>_gaugetensor.csv (nsp)
    X = table();
    for k = 1:numel(f)
        if f(k).name == "gaugetensor.csv", [~, sub] = fileparts(f(k).folder); pre = ''; else
            sub = extractBefore(f(k).name, "_gaugetensor.csv"); pre = [sub '_']; end
        t = readtable(fullfile(f(k).folder, f(k).name), 'TextType', 'string');  t.hemi = string(t.hemi);
        ds = "";  m = fullfile(f(k).folder, [pre 'metrics.csv']);
        if isfile(m)
            mm = readtable(m, 'TextType', 'string', 'Delimiter', ',');
            if ismember('dataset', mm.Properties.VariableNames) && height(mm), ds = string(mm.dataset(1)); end
        end
        t.subject = repmat(string(sub), height(t), 1);  t.dataset = repmat(ds, height(t), 1);
        X = [X; t]; %#ok<AGROW>
    end
    if isempty(X), error('scale:reducegauge:empty', 'no gaugetensor.csv under %s', inDir); end
    X.dataset(ismissing(X.dataset)) = "";
    [g, G] = findgroups(X(:, {'dataset','hemi','depth','node_id'}));
    G.n = splitapply(@numel, X.t11, g);
    G.t11 = splitapply(@(x) mean(x, 'omitnan'), X.t11, g);
    G.t22 = splitapply(@(x) mean(x, 'omitnan'), X.t22, g);
    G.t12 = splitapply(@(x) mean(x, 'omitnan'), X.t12, g);
    G.anisotropy = sqrt((G.t11 - G.t22).^2 + 4*G.t12.^2) ./ (G.t11 + G.t22);
    G.axis_deg = 0.5 * atan2d(2*G.t12, G.t11 - G.t22);
    G.axis_R = splitapply(@(a) abs(mean(exp(2i*deg2rad(a)), 'omitnan')), X.axis_deg, g);
    G.anisotropy_subject_median = splitapply(@(x) median(x, 'omitnan'), X.anisotropy, g);
    G.normal_share_median = splitapply(@(x) median(x, 'omitnan'), X.normal_share, g);
    G.colat_median_deg = splitapply(@(x) median(x, 'omitnan'), X.colat_median_deg, g);
    G.polar = G.colat_median_deg < opts.PolarDeg | G.colat_median_deg > 180 - opts.PolarDeg;
    G.axis_deg(G.polar) = NaN;  G.axis_R(G.polar) = NaN;
    writetable(G, fullfile(outDir, 'gaugegroup.csv'));
    fprintf('[reducegauge] %d subjects, %d tile rows (%d polar)\n', numel(unique(X.subject)), height(G), nnz(G.polar));
end

% Author: Diellor Basha, 2026
