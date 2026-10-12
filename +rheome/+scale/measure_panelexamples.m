function [T, Ex] = measure_panelexamples(name, S)
% SCALE.MEASURE_PANELEXAMPLES  The single-participant example data of MS1 Figs. 4-7 panel A (figure data, no group metric).
%
%   [T, Ex] = rheome.scale.measure_panelexamples(name)
%
% One example per figure, each made by the measure whose group numbers the figure shows, with small options
% (request 64b475b1, for Graphics 059029b9):
%   fig4  rheome.scale.measure_plantfloors, one source + one vortex at 130 mm on the left hemisphere, noise-free,
%         10 dB and without the instrument: the planted current and the recovered band maps
%   fig5  the point spread of the whitened plain MNE for the left-hemisphere seed whose r50 is the median of the
%         300 seeds rheome.scale.measure_resolution draws (rheome.inverse.resolution, seed 7), with every seed's r50
%   fig6  rheome.scale.measure_inject, one 0.05 m/s mover at depth 7 and 18.75 fps (JLevel 4): the planted path
%         and the tracked one
%   fig7  rheome.scale.measure_aperiodic, the rest arm, an offset x100 patch (sigma 100 mm) moving at 0.15 m/s:
%         the patch on the cortex and the sensor spectra (background, plant, total, periodic split)
% The Fig. 2 example (per-window path scores) is the grouptrack analysis's grouptrackwindows.csv.
% T: rheome.scale.rows (analysis "panelexamples"): each example's headline readout, for a check against the group.
%
% ⚠ Examples, not estimates: one plant each, drawn with the measure's own seed at these options -- so not the
%   same draw as the group run's first plant.
%
% See also: rheome.scale.run, rheome.scale.measure_plantfloors, rheome.inverse.resolution,
%           rheome.scale.measure_inject, rheome.scale.measure_aperiodic
%
% Author: Diellor Basha, 2026

    if nargin < 2 || isempty(S), S = rheome.scale.sensors(name); end
    Ex = struct();
    [~, ~, Ex.fig4] = rheome.scale.measure_plantfloors(name, S, LambdaMM=130, PlantsPerCell=1, Hemis="L", SNRdB=[Inf 10]);
    Ex.fig5 = i_psf(S);
    [~, ~, Ex.fig6] = rheome.scale.measure_inject(name, S, Depths=7, JLevels=4, Cases="move05", Speeds=0.05, NPer=2, Hemis="L");
    [~, ~, Ex.fig7] = rheome.scale.measure_aperiodic(name, S, Arms="rest", SigmaMM=100, OffsetX=100, SpeedsMS=0.15, ...
                                                     FitWindowsS=[6 Inf], Depths=7, Reps=1, NullWindows=8);
    Ex = i_plain(Ex);                    % char and structs only: the -v7 file must load in scipy.io.loadmat
    m = ["fig4_source_err_mm_10dB" "fig4_vortex_err_mm_10dB" "fig5_r50_mm" "fig5_r50_median_mm" "fig6_found" ...
         "fig6_speed_ms" "fig7_plant_inband_dB"];
    k = find(Ex.fig4.snr_dB == 10, 1);
    v = [Ex.fig4.source_err_mm(k) Ex.fig4.vortex_err_mm(k) Ex.fig5.r50_mm Ex.fig5.r50_median_mm ...
         Ex.fig6.readout.found Ex.fig6.readout.speedMS Ex.fig7.plant_inband_dB];
    T = rheome.scale.rows("panelexamples", m, v, ["mm" "mm" "mm" "mm" "flag" "m/s" "dB"]);
end

function x = i_plain(x)
% strings -> char (cellstr if not scalar), tables -> structs of columns, recursively
    if isstring(x), if isscalar(x), x = char(x); else, x = cellstr(x); end
    elseif istable(x), x = i_plain(table2struct(x, 'ToScalar', true));
    elseif isstruct(x)
        for i = 1:numel(x), for f = string(fieldnames(x))', x(i).(f) = i_plain(x(i).(f)); end, end
    end
end

function P = i_psf(S)
% the orientation-free PSF power of one seed, as rheome.inverse.resolution's localisation family computes it
    H = S.B.L;  SL = H.S;  K = S.Res.ImagingKernel;
    out = rheome.inverse.resolution(K, S.G, H.lbo, Vertices=SL.Vertices, Faces=SL.Faces, ...
              Normals=SL.VertNormals, GlobalIdx=H.gv, SensorLoc=S.Loc, NumSeeds=300);
    [~, j] = min(abs(out.r50 - median(out.r50)));  sd = out.seeds(j);  gv = double(H.gv(:));
    Eb = K * S.G(:, (gv(sd)-1)*3 + (1:3));  p = sum(reshape(sum(Eb.^2, 2), 3, []), 1)';  pL = p(gv);
    V = SL.Vertices;  F = double(SL.Faces);  E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    d = distances(simplify(graph(E(:,1), E(:,2), vecnorm(V(E(:,1),:) - V(E(:,2),:), 2, 2), size(V,1))), sd)';
    q = pL / sum(pL);  [ds, o] = sort(d);  r50 = ds(find(cumsum(q(o)) >= 0.5, 1));  [~, im] = max(pL);
    assert(abs(r50 - out.r50(j)) < 1e-9, 'scale:panelexamples:psf', 'the PSF does not reproduce rheome.inverse.resolution');
    P = struct('hemi', "L", 'vertices_mm', 1e3*V, 'faces', F, 'source_vertex', sd, 'psf_power', pL, ...
        'geodesic_mm', 1e3*d, 'r50_mm', 1e3*r50, 'ple_mm', 1e3*d(im), 'peak_vertex', im, ...
        'offhemi_fraction', 1 - sum(pL)/sum(p), 'r50_median_mm', 1e3*median(out.r50), 'lambda_loc_mm', 1e3*out.lambdaLoc, ...
        'seed_vertex', out.seeds(:)', 'seed_r50_mm', 1e3*out.r50(:)', 'seed_ple_mm', 1e3*out.ple(:)', ...
        'seed_depth_mm', 1e3*out.depth(:)', 'inverse', "whitened plain MNE, SnrFixed 3 (rheome.scale.sensors)");
end

% Author: Diellor Basha, 2026
