function p = window(p, span, opts)
% PIPELINE.WINDOW  Append a selection on the frame axis: the bookkeeping step.
%
%   p = rheome.pipeline.window(p, [262 270], Rate=600)        % seconds
%   p = rheome.pipeline.window(p, [1000 5800], Units="samples")
%
% ⭐ A SELECTION IS A STEP, NOT A PREAMBLE. Which stretch of the recording an analysis ran
% on is part of what produced the answer, so it belongs in the plan beside the operators and
% travels into the record. The field type does not change; the frame domain does, and the
% new one remembers the span it came from -- which is the same key the tile store labels and
% measurements use (level, k), so a window chosen by a query and a window chosen by hand are
% recorded the same way.
%
% Author: Diellor Basha, 2026

    arguments
        p (1,1) struct
        span (1,2) double
        opts.Rate  (1,1) double = NaN
        opts.Units (1,1) string {mustBeMember(opts.Units, ["seconds","samples"])} = "seconds"
        opts.Note  (1,1) string = ""
    end
    fs = opts.Rate;
    if opts.Units == "seconds" && isnan(fs)
        error('pipeline:window:rate', 'A window in seconds needs Rate (the sampling rate).');
    end
    if opts.Units == "seconds"
        a = max(1, floor(span(1)*fs) + 1);  b = ceil(span(2)*fs);
    else
        a = max(1, round(span(1)));  b = round(span(2));
    end
    if b < a, error('pipeline:window:span', 'An empty window: %g to %g.', span(1), span(2)); end
    n = b - a + 1;
    fr = rheome.domain.index(n, Name=sprintf('%s[%d:%d]', p.frames.name, a, b), Units=p.frames.units);
    s = struct('kind', 'window', 'id', "window", 'in', p.type, 'out', p.type, ...
               'domain_in', p.domain, 'domain_out', p.domain, 'frames_out', fr, 'range', [a b], ...
               'note', sprintf('samples %d:%d (%g %s)%s', a, b, diff(span), opts.Units, i_note(opts.Note)));
    p.frames = fr;
    p.steps(end+1) = s;
end

function s = i_note(n)
    if n == "", s = ''; else, s = sprintf(' -- %s', char(n)); end
end
% Author: Diellor Basha, 2026
