function [R, costs] = cohort(C, fcn)
% SELECT.COHORT  Run a query over every store of a catalogue selection and stack the results.
%
%   [R, costs] = rheome.select.cohort(C, @(db) rheome.select.frames(db, q))
%
% Relative thresholds resolve per recording inside fcn (the spec's rule: absolute MEG
% amplitudes vary with head position). costs is a struct array, one per store.
%
% Author: Diellor Basha, 2026

    parts = cell(height(C), 1);  costs = cell(height(C), 1);
    for i = 1:height(C)
        db = rheome.select.open(char(C.file(i)));
        [Ri, ci] = fcn(db);
        if ~ismember('recording_id', Ri.Properties.VariableNames), Ri.recording_id = repmat(string(db.recording_id), height(Ri), 1); end
        parts{i} = Ri;  costs{i} = ci;
    end
    R = vertcat(parts{:});
    costs = [costs{:}];
end
% Author: Diellor Basha, 2026
