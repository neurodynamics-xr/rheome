function fmri = bold(files, nV)
% IO.READ.BOLD  Resting-state BOLD on the cortex from nsp bst-fmri's Brainstorm overlay, one entry per run.
%
%   fmri = rheome.io.read.bold(files)        % files: cell of .../fmri_<ses>_task-rest/results_surface_fmri_*.mat
%   fmri = rheome.io.read.bold(files, nV)    % and check every map has nV vertices
%
% Each results file (ImageGridAmp [nV x nT], Time, SurfaceFile, DisplayUnits) is read with its run's
% fMRIPrep confounds, matrix_confounds_fmri_<tag>.mat in the same folder (Value [nCol x nT], Description).
% Returns a struct array: .name (task_rest_run_1 ...) .x [nV x nT] single .tr (s) .units .confounds [nT x nCol]
% .confoundNames (string) .file.
%
% ⭐ The values are fMRIPrep's preprocessed BOLD (a.u.) on fsaverage5, interpolated to the subject's cortex
% (Shepard, 8 neighbours, registered spheres): no filtering, no confound regression, no smoothing. The cortex
% therefore carries no detail finer than fsaverage5's spacing (~10k vertices per hemisphere).
% ⚠ A map on another surface than tess_cortex_pial_low (the leadfield's) is an error: it would mis-index
% every vertex.
%
% See also: rheome.scale.importsubject, rheome.scale.measure_multimodal
%
% Author: Diellor Basha, 2026

    if nargin < 2, nV = []; end
    fmri = struct('name', {}, 'x', {}, 'tr', {}, 'units', {}, 'confounds', {}, 'confoundNames', {}, 'file', {});
    for k = 1:numel(files)
        f = files{k};
        r = builtin('load', f, 'ImageGridAmp', 'Time', 'SurfaceFile', 'DisplayUnits');
        [~, sf] = fileparts(char(r.SurfaceFile));
        if ~strcmp(sf, 'tess_cortex_pial_low')
            error('io:read:bold:surface', '%s is on %s, not tess_cortex_pial_low', f, r.SurfaceFile);
        end
        if ~isempty(nV) && size(r.ImageGridAmp, 1) ~= nV
            error('io:read:bold:vertices', '%s has %d vertices, the cortex %d', f, size(r.ImageGridAmp, 1), nV);
        end
        nT = size(r.ImageGridAmp, 2);  dt = diff(double(r.Time));
        if nT < 2 || numel(r.Time) ~= nT || any(abs(dt - dt(1)) > 1e-6 * dt(1)) || dt(1) <= 0
            error('io:read:bold:time', '%s: Time is not a regular axis of %d volumes', f, nT);
        end
        if any(~isfinite(r.ImageGridAmp(:))), error('io:read:bold:finite', '%s has non-finite values', f); end
        [d, b] = fileparts(f);  tag = erase(b, 'results_surface_fmri_');
        cf = fullfile(d, ['matrix_confounds_fmri_' tag '.mat']);
        if ~isfile(cf), error('io:read:bold:confounds', 'no %s beside %s', cf, f); end
        c = builtin('load', cf, 'Value', 'Description');
        if size(c.Value, 2) ~= nT, error('io:read:bold:confounds', '%s has %d samples, the BOLD %d', cf, size(c.Value, 2), nT); end
        fmri(end+1) = struct('name', tag, 'x', single(r.ImageGridAmp), 'tr', dt(1), 'units', char(r.DisplayUnits), ...
                             'confounds', double(c.Value)', 'confoundNames', string(c.Description(:))', 'file', f); %#ok<AGROW>
    end
end

% Author: Diellor Basha, 2026
