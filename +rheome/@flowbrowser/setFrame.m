function setFrame(app, i)
% SETFRAME  Move to a frame: clamp, sync the slider and the readout, refresh.
%
% The single entry point for changing time, so the slider, the step buttons, the keyboard
% and the play loop cannot drift out of sync with each other.
%
% Author: Diellor Basha, 2026

    i = min(max(round(i), 1), app.Page.NumFrames);
    if i == app.Frame && ~isempty(app.FrameLbl.Text), return; end
    app.Frame = i;
    if isvalid(app.FrameSlider), app.FrameSlider.Value = i; end   % programmatic: no callback
    app.updateFrameLabel();
    app.refresh();
end

% Author: Diellor Basha, 2026
