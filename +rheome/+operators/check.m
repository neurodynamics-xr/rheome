function [out, op] = check(id, data, d, opts)
% OPERATORS.CHECK  Is this field admissible as that operator's input, and what comes back?
%
%   [outType, op] = rheome.operators.check("weak_curl", J, dCortex)
%   outType = rheome.operators.check("weak_curl", [], dCortex)        % just ask for the arrow
%
% ⭐ THE QUESTION AN ANALYSIS KEEPS ASKING. Given a field and an operator, does it fit, and
% what am I holding afterwards? The answer is the registry's arrow, checked against the
% field's scalar count on this domain, so a chain can be assembled and verified before any
% of it runs.
%
% Returns the OUTPUT field type id, and the operator's registry row.
%
% Author: Diellor Basha, 2026

    arguments
        id
        data = []
        d (1,1) struct = struct('kind','complex','nV',0,'nE',0,'nF',0)
        opts.Throw (1,1) logical = true
    end
    T = rheome.operators.registry();
    i = find(T.id == string(id), 1);
    if isempty(i)
        error('operators:check:unknown', 'No operator ''%s'' in the registry.', char(string(id)));
    end
    op = table2struct(T(i, :));
    if op.status == "planned" && opts.Throw
        error('operators:check:planned', ...
              '%s is declared but not built here (%s). rheome.fieldtype.registry marks the gap.', op.id, op.fcn);
    end
    if ~isempty(data)
        rheome.fieldtype.validate(data, op.in, d, Throw=opts.Throw, Name=sprintf('input of %s', op.id));
    end
    out = op.out;
end
% Author: Diellor Basha, 2026
