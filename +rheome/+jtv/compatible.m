function [ok, why] = compatible(a, b, opts)
% JTV.COMPATIBLE  Do two joint representations live on the same axes?
%
%   [ok, why] = rheome.jtv.compatible(axA, axB)
%
% Checks the eigenvalue axis, the frequency axis, and the conventions that make two
% representations combinable at all. Call it before composing, comparing or reconstructing across
% two analyses -- particularly across BANDS, since a per-band run truncates to its own bins and two
% such runs do NOT share a frequency axis even though both are called "the joint spectrum".
%
% ⚠ .half and .boundary are checked because a mismatch there is silent and fatal. A positive-half
% (analytic, complex) representation and a full-axis (real) one cannot be combined; nor can a
% ring-boundary and a path-boundary one, whose implied extensions differ.
%
% INPUTS:  a, b  axes structs (rheome.jtv.axes)   opts .tol relative tolerance (default 1e-10)
% OUTPUT:  ok    logical                   why  cellstr of the mismatches found
%
% See also: rheome.jtv.axes
%
% Author: Diellor Basha, 2026

    if nargin < 3, opts = struct(); end
    if ~isfield(opts,'tol') || isempty(opts.tol), opts.tol = 1e-10; end
    why = {};

    if isfield(a,'lambda') && isfield(b,'lambda')
        if numel(a.lambda) ~= numel(b.lambda)
            why{end+1} = sprintf('lambda: %d vs %d modes', numel(a.lambda), numel(b.lambda));
        elseif max(abs(a.lambda - b.lambda)) > opts.tol * max(a.lmax, eps)
            why{end+1} = 'lambda: same count, different values (different anatomy or basis)';
        end
    end
    if isfield(a,'f') && isfield(b,'f')
        if numel(a.f) ~= numel(b.f)
            why{end+1} = sprintf('frequency: %d vs %d bins -- different bands retained?', ...
                numel(a.f), numel(b.f));
        elseif max(abs(a.f - b.f)) > opts.tol * max(max(a.f), eps)
            why{end+1} = 'frequency: same count, different values';
        end
    end
    if isfield(a,'half') && isfield(b,'half') && ~strcmpi(a.half, b.half)
        why{end+1} = sprintf('half: ''%s'' vs ''%s'' -- one is analytic/complex, the other real', ...
            a.half, b.half);
    end
    if isfield(a,'boundary') && isfield(b,'boundary') && ~strcmpi(a.boundary, b.boundary)
        why{end+1} = sprintf('boundary: ''%s'' vs ''%s'' -- different implied time extension', ...
            a.boundary, b.boundary);
    end
    ok = isempty(why);
end

% Author: Diellor Basha, 2026
