function updateFrameLabel(app)
% UPDATEFRAMELABEL  Frame number, time, and whether we are in the page's margin.
%
% ⚠ The margin is called out here rather than only shaded on the Timeline, because the
% Cortex tab has no time axis to shade -- without this you cannot tell from the surfaces
% alone that you are looking at cone-contaminated frames.
%
% Author: Diellor Basha, 2026

    fp = app.Page;
    tag = '';
    if ~fp.Core(app.Frame), tag = '  ⚠ margin'; end
    app.FrameLbl.Text = sprintf('%d/%d   %.3f s%s', ...
        app.Frame, fp.NumFrames, fp.Time(app.Frame), tag);
end

% Author: Diellor Basha, 2026
