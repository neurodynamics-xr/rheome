function [Y, rec] = run(p, e, X, opts)
% PIPELINE.RUN  Execute a plan and return the result WITH its bookkeeping.
%
%   [g, rec] = rheome.pipeline.run(p, e, F);
%
% ⭐ THE RECORD IS THE POINT. What came out is a matrix like any other; what makes it usable
% a week later is knowing which recording and window it came from, which operators were
% applied in which order, which domain each intermediate lived on, and what shape it had.
% `rec.steps` is one row per step with all of that, and `rec.out` is the field type and
% domain of the answer -- so the result can be validated, stored, or handed to the next
% pipeline without guessing.
%
% ⚠ EVERY STEP IS VALIDATED AS IT RUNS, not only at the end: the intermediate is checked
% against the type the registry promised on the domain the environment resolved. A wrong
% applier shows up at the step that broke the contract, not as a mysterious final shape.
%
% Author: Diellor Basha, 2026

    arguments
        p (1,1) struct
        e (1,1) struct
        X
        opts.Validate (1,1) logical = true
    end
    T = rheome.operators.registry();
    rheome.fieldtype.validate(X, p.in.type, p.in.domain, Name="pipeline input");
    d = p.in.domain;  fr = p.in.frames;  t = p.in.type;
    Y = X;
    rec = struct('name', p.name, 'env', e.name, 'in', p.in, 'started', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
                 'steps', struct('step', {}, 'kind', {}, 'id', {}, 'in', {}, 'out', {}, ...
                                 'domain', {}, 'frames', {}, 'rows', {}, 'cols', {}, 'seconds', {}, 'note', {}));
    for i = 1:numel(p.steps)
        s = p.steps(i);  t0 = tic;
        if s.kind == "window"
            b = min(s.range(2), size(Y, 2));
            Y = Y(:, s.range(1):b);
            fr = s.frames_out;
        else
            op = table2struct(T(T.id == string(s.id), :));
            fn = i_fn(op.apply);          % MATLAB cannot call a call's result
            Y = fn(e, Y);
            t = op.out;
            if op.target ~= "", d = e.domains.(char(op.target)); end
            if opts.Validate
                rheome.fieldtype.validate(Y, t, d, Name=sprintf('output of step %d (%s)', i, s.id));
            end
        end
        rec.steps(end+1) = struct('step', i, 'kind', s.kind, 'id', string(s.id), 'in', string(s.in), ...
            'out', string(t), 'domain', string(i_dname(d)), 'frames', string(fr.name), ...
            'rows', size(Y, 1), 'cols', size(Y, 2), 'seconds', toc(t0), 'note', string(s.note));
    end
    rec.out = struct('type', string(t), 'domain', d, 'frames', fr);
    rec.table = struct2table(rec.steps, 'AsArray', true);
end

function s = i_dname(d)
    if isempty(d.name), s = sprintf('%s/%s', d.kind, d.id); else, s = d.name; end
end

% table2struct hands a cell column back as a cell for some rows and bare for others; one
% unwrap here keeps both call sites from caring.
function f = i_fn(a)
    if iscell(a), f = a{1}; else, f = a; end
end
% Author: Diellor Basha, 2026
