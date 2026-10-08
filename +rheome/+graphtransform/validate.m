function [ok, why] = validate(T)
% GRAPHTRANSFORM.VALIDATE  Is T a usable graphfilterbank transform?
%
%   [ok, why] = rheome.graphtransform.validate(T)
%
% See also: rheome.graphtransform.eigen, rheome.graphtransform.chebyshev
%
% Author: Diellor Basha, 2026

    why = {};
    if ~isstruct(T)
        why{end+1} = 'not a struct';  ok = false;  return;
    end
    for f = {'forward','inverse','norm'}
        if ~isfield(T, f{1}) || ~isa(T.(f{1}), 'function_handle')
            why{end+1} = sprintf('missing or non-handle .%s', f{1});
        end
    end
    if isempty(why) && isfield(T, 'lambda') && ~isempty(T.lambda) && isfield(T, 'rows')
        F = randn(T.rows, 2);
        if ~isequal(size(T.forward(F)), [numel(T.lambda), 2])
            why{end+1} = 'forward(F) does not return [numel(lambda) x nT]';
        end
    end
    ok = isempty(why);
end

% Author: Diellor Basha, 2026
