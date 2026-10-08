function T = describe(p)
% PIPELINE.DESCRIBE  The plan as a table: every step, its arrow, and where it lands.
%
%   disp(rheome.pipeline.describe(p))
%
% Needs no environment and runs nothing: this is what a plan IS. The same table is what a
% reviewer reads and what a record is compared against afterwards.
%
% Author: Diellor Basha, 2026

    n = numel(p.steps);
    step = (1:n)';  kind = strings(n,1);  id = strings(n,1);
    in = strings(n,1);  out = strings(n,1);  dom = strings(n,1);  fr = strings(n,1);  note = strings(n,1);
    for i = 1:n
        s = p.steps(i);
        kind(i) = string(s.kind);  id(i) = string(s.id);
        in(i) = string(s.in);  out(i) = string(s.out);
        if isstruct(s.domain_out), dom(i) = string(i_dname(s.domain_out)); else, dom(i) = string(s.domain_out); end
        fr(i) = string(s.frames_out.name);  note(i) = string(s.note);
    end
    T = table(step, kind, id, in, out, dom, fr, note, ...
              'VariableNames', {'step','kind','operator','in','out','domain','frames','note'});
end

function s = i_dname(d)
    if isempty(d.name), s = sprintf('%s/%s', d.kind, d.id); else, s = d.name; end
end
% Author: Diellor Basha, 2026
