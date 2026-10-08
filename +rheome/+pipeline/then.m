function p = then(p, id)
% PIPELINE.THEN  Append an operator, checking the arrow against the state right now.
%
%   p = rheome.pipeline.then(p, "inverse_mne");
%   p = rheome.pipeline.then(p, "curl");
%
% ⭐ THE CHECK IS AT BUILD TIME. The registry says what this operator consumes; the pipeline
% knows what it is holding. A mismatch is refused HERE, naming both types, rather than
% surfacing as a shape error inside a matrix product three steps later.
%
% ⚠ THE DOMAIN MOVES WITH THE ARROW. A cross-domain operator names which of the
% environment's domains its output lives on (`target`), so after inverse_mne the pipeline is
% on the cortex and after lb_forward it is on the mode line -- and the operators offered
% next are the ones that domain supports.
%
% Author: Diellor Basha, 2026

    T = rheome.operators.registry();
    i = find(T.id == string(id), 1);
    if isempty(i)
        error('pipeline:then:unknown', 'No operator ''%s'' in rheome.operators.registry.', char(string(id)));
    end
    op = table2struct(T(i, :));
    if op.status ~= "built"
        error('pipeline:then:planned', '%s is declared but not built here (%s).', op.id, op.fcn);
    end
    if op.in ~= p.type
        error('pipeline:then:arrow', ...
              ['%s takes a %s; the pipeline is holding a %s. Operators that accept it: %s.'], ...
              op.id, op.in, p.type, strjoin(cellstr(T.id(T.in == p.type & T.status == "built")), ', '));
    end
    s = struct('kind', 'operator', 'id', op.id, 'in', op.in, 'out', op.out, ...
               'domain_in', p.domain, 'domain_out', [], 'frames_out', p.frames, ...
               'range', [], 'note', char(op.fcn));   % range: the field a window step fills
    p.type = op.out;
    if op.target ~= ""                                  % the domain changes; run() resolves it
        s.domain_out = char(op.target);
        p.domain = struct('kind', 'pending', 'nV', NaN, 'nE', NaN, 'nF', NaN, 'id', '', ...
                          'name', char(op.target), 'source', '');
    else
        s.domain_out = p.domain;
    end
    p.steps(end+1) = s;
end
% Author: Diellor Basha, 2026
