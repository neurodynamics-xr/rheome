function c = fbd_check(name, measured, expected, tol, unit)
% FBD_CHECK  One self-check: compare a measured quantity against its analytic value.
%
%   c = fbd_check('bump 1 width', 20.4, 20.0, 1.6, 'mm')
%
% Prints one aligned PASS/FAIL line and returns the record, so a demo both reports as it
% runs and can be asserted on afterwards.
%
% ⚠ A NaN measurement is a FAILURE. Returning NaN from a fit that could not converge and
% letting abs(NaN - x) <= tol evaluate to false is right, but only by accident -- this
% states it.
%
% Author: Diellor Basha, 2026

    if nargin < 5, unit = ''; end
    ok = isfinite(measured) && isfinite(expected) && abs(measured - expected) <= tol;
    c = struct('name', string(name), 'measured', measured, 'expected', expected, ...
               'tol', tol, 'unit', string(unit), 'pass', ok);
    fprintf('  %-38s %10.4g vs %-10.4g (+-%.3g) %-6s %s\n', ...
        name, measured, expected, tol, unit, fbd_ternary(ok, 'PASS', 'FAIL'));
end

% Author: Diellor Basha, 2026
