function [f, why] = compile(p, e, opts)
% PIPELINE.COMPILE  Can this plan run in this environment, and as what function?
%
%   f = rheome.pipeline.compile(p, e);   Y = f(X);
%   [f, why] = rheome.pipeline.compile(p, e, Throw=false)
%
% ⭐ EVERY REASON IT COULD FAIL, BEFORE IT RUNS. A step whose operator has no applier yet, a
% cross-domain step whose target domain the environment does not hold, an environment
% missing the matrix an applier reaches for: all of them are found here, named, and reported
% together rather than one per run.
%
% The returned handle takes the input field and gives the output field. rheome.pipeline.run wraps
% it and also returns the record.
%
% Author: Diellor Basha, 2026

    arguments
        p (1,1) struct
        e (1,1) struct
        opts.Throw (1,1) logical = true
    end
    why = string.empty(0, 1);
    T = rheome.operators.registry();
    need = struct('mass', 'M', 'laplace_beltrami', 'L', 'inverse_mne', 'K');
    for i = 1:numel(p.steps)
        s = p.steps(i);
        if s.kind ~= "operator", continue; end
        op = table2struct(T(T.id == string(s.id), :));
        if isempty(i_fn(op.apply))
            why(end+1) = sprintf('step %d (%s): no applier yet; %s builds it but nothing wires it', ...
                                 i, s.id, op.fcn);                                  %#ok<AGROW>
        end
        if op.target ~= "" && ~isfield(e.domains, char(op.target))
            why(end+1) = sprintf('step %d (%s): the environment has no ''%s'' domain', ...
                                 i, s.id, op.target);                               %#ok<AGROW>
        end
        for fn = string(fieldnames(need))'
            if s.id == fn && isempty(e.(need.(char(fn))))
                why(end+1) = sprintf('step %d (%s): the environment has no %s', i, s.id, need.(char(fn)));  %#ok<AGROW>
            end
        end
        if any(strcmp(s.id, {'curl','divergence','weak_curl','weak_div','face_gradient','face_average','poisson'})) ...
                && (isempty(e.fg) || isempty(e.S))
            why(end+1) = sprintf('step %d (%s): the environment has no surface operators', i, s.id);  %#ok<AGROW>
        end
        if any(strcmp(s.id, {'lb_forward','lb_inverse'})) && isempty(e.lbo)
            why(end+1) = sprintf('step %d (%s): the environment has no eigenbasis', i, s.id);  %#ok<AGROW>
        end
    end
    if ~isempty(why)
        if opts.Throw
            error('pipeline:compile:unbound', 'This plan cannot run here:\n  %s', strjoin(cellstr(why), sprintf('\n  ')));
        end
        f = [];  return
    end
    f = @(X) i_apply(p, e, X);
end

function Y = i_apply(p, e, X)
    T = rheome.operators.registry();
    Y = X;
    for i = 1:numel(p.steps)
        s = p.steps(i);
        if s.kind == "window"
            Y = Y(:, s.range(1):min(s.range(2), size(Y, 2)));
        else
            op = table2struct(T(T.id == string(s.id), :));
            fn = i_fn(op.apply);          % MATLAB cannot call a call's result
            Y = fn(e, Y);
        end
    end
end

% table2struct hands a cell column back as a cell for some rows and bare for others; one
% unwrap here keeps both call sites from caring.
function f = i_fn(a)
    if iscell(a), f = a{1}; else, f = a; end
end
% Author: Diellor Basha, 2026
