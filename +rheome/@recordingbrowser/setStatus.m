function setStatus(app, s)
% SETSTATUS  One line of status beside the readout.
% Author: Diellor Basha, 2026
    if ~isempty(app.StatusLbl) && isvalid(app.StatusLbl), set(app.StatusLbl, 'String', s); end
end
% Author: Diellor Basha, 2026
