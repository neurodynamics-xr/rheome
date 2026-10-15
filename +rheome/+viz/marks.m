function file = marks(m, render)
% VIZ.MARKS  Write a render's marks sidecar (<render>.marks.json): every overlay drawn on the cortex, with its key.
%
%   file = rheome.viz.marks(out.marks, 'figs/F2b.png')    % -> figs/F2b.marks.json
%
% The viewer's screenshot and figure events carry 'marks' (format nxr.marks@1: the surface, each overlay's style and
% key); rheome.view's screenshot and rheome.snapshot write it beside each PNG with this. It is what figure_check.py E1
% reads (Chekhov's gun, the figure tokens' cortex.marks): an overlay without a key fails the figure.
%
% INPUTS:
%   m       the marks struct (decoded JSON) with .format 'nxr.marks@1'
%   render  the render's file (its extension is replaced) or its path without one
%
% Example:
%   f = rheome.viz.marks(struct('format', 'nxr.marks@1', 'overlays', {{}}), fullfile(tempdir, 'r.png'));
%   assert(endsWith(f, 'r.marks.json')); delete(f)
%
% See also: rheome.view, rheome.snapshot
%
% Author: Diellor Basha, 2026
    if ~isstruct(m) || ~isfield(m, 'format') || ~strcmp(m.format, 'nxr.marks@1')
        error('rheome:viz:marks', 'Marks of format nxr.marks@1 are expected (the viewer''s screenshot or figure event).');
    end
    [d, n] = fileparts(char(render));
    file = fullfile(d, [n '.marks.json']);
    fid = fopen(file, 'w');
    if fid < 0, error('rheome:viz:marks', 'Cannot write %s.', file); end
    fwrite(fid, jsonencode(m, 'PrettyPrint', true));  fclose(fid);
end
