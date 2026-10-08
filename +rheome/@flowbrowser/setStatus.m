function setStatus(app, s)
% SETSTATUS  One-line status in the control strip.
% Author: Diellor Basha, 2026
    if ~isempty(app.StatusLbl) && isvalid(app.StatusLbl)
        app.StatusLbl.Text = s;
    end
end

% Author: Diellor Basha, 2026
