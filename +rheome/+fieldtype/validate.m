function [ok, why] = validate(data, id, d, opts)
% FIELDTYPE.VALIDATE  Does this array hold a field of this type on this domain?
%
%   rheome.fieldtype.validate(J, "ambientVertexWorld", dCortex)          % throws on a mismatch
%   [ok, why] = rheome.fieldtype.validate(J, "ambientVertexWorld", dCortex, Throw=false)
%
% ⭐ THE CHECK IS A SCALAR COUNT, NOT A ROW COUNT. A field of type t on domain d has
% elements(d, t.rank) * components(t) scalars per frame -- an ambient vertex field on 20484
% vertices is 61452 rows, not 20484, and an [nV x 3] face gradient is 3*nF scalars laid out
% the other way round. Checking rows alone is what lets a constrained leadfield pass as an
% unconstrained one; rheome.inverse.mne documents that exact collision, which is why it has to ask
% for the vertex count separately.
%
% ⚠ THE TRAILING DIMENSION IS FREE BY DEFAULT. Columns are frames (time, trials), which are a separate
% domain: the type says nothing about how many there are. ⭐ Pass Frames=n -- or a product domain's
% .nFrames -- to check it too. That is opt-in because the default looseness is what every existing
% caller relies on, and because turning it on everywhere at once surfaces every sloppy frame count in
% one pass rather than where someone is ready to fix it.
%
% Author: Diellor Basha, 2026

    arguments
        data
        id
        d (1,1) struct
        opts.Throw  (1,1) logical = true
        opts.Name   (1,1) string = "field"
        opts.Frames        double = []
    end
    why = '';
    f = rheome.fieldtype.describe(id);
    if ~ismember(string(d.kind), string(f.domain_kind))
        why = sprintf('%s is a %s field; the domain is a %s', f.id, strjoin(cellstr(f.domain_kind), '/'), d.kind);
    else
        n = rheome.domain.elements(d, f.rank);
        want = n * f.components;
        got = size(data, 1);
        if got ~= want && ~(f.components > 1 && isequal(size(data), [n f.components]))
            why = sprintf(['%s on this domain needs %d scalars per frame (%d %s elements x %d ' ...
                           'components), got %d rows'], f.id, want, n, f.rank, f.components, got);
        end
    end
    % ⭐ OPT-IN FRAME CHECK. The trailing dimension is free by default and always has been, which is
    % why a wrong frame count is invisible to the type system -- rheome.filters.impulse in 'js' mode returned
    % the wrong number of samples in rheome.forward.simulate and nothing noticed until the arithmetic failed
    % on `m*B1 + W{k}`. Pass Frames (or a product domain's .nFrames) to close that hole at the type.
    % ⚠ DELIBERATELY OPT-IN: turning it on everywhere at once would surface every existing looseness
    % about frame counts in one go, so callers adopt it where they mean it.
    if isempty(why) && ~isempty(opts.Frames)
        gotF = size(data, 2);
        if gotF ~= opts.Frames
            why = sprintf('%s should hold %d frames, got %d columns', f.id, opts.Frames, gotF);
        end
    end
    ok = isempty(why);
    if ~ok && opts.Throw
        error('fieldtype:validate:shape', '%s: %s.', char(opts.Name), why);
    end
end
% Author: Diellor Basha, 2026
