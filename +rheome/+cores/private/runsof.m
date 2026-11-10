function r = runsof(mask)
% RUNSOF  Runs of true in a logical vector, as 0-based [start stop) pairs, one row per run (nxr runs_of).
    d = diff([0, double(mask(:).'), 0]);
    r = [find(d == 1).' - 1, find(d == -1).' - 1];
end
