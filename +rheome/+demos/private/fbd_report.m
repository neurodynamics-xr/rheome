function [ok, n] = fbd_report(checks)
% FBD_REPORT  Summarise a run of checks.
%   [ok, n] = fbd_report(checks)
% Author: Diellor Basha, 2026
    n  = numel(checks);
    np = sum([checks.pass]);
    ok = (np == n);
    fprintf('  ---------------------------------------------------------------\n');
    fprintf('  %d/%d checks %s\n', np, n, fbd_ternary(ok, 'PASS', 'FAIL'));
end

% Author: Diellor Basha, 2026
