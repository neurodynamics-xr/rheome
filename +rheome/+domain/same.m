function [ok, why] = same(dA, dB)
% DOMAIN.SAME  Are two fields on the same domain? Ids for simple domains, factors and conventions
% for a product.
%
%   [ok, why] = rheome.domain.same(dA, dB)
%
% ⭐ THE ONE EQUALITY. A domain "carries no coordinates: two fields are on the same domain when the
% ids match" (+domain/Contents.m), and for a product that has to extend to the factors AND the
% CONVENTIONS -- two banks can agree on every count and be incompatible because one holds the
% positive half of the spectrum and the other the whole of it, or one treats time as a ring and the
% other as a path. @jointfilterbank/iscompatible existed to catch exactly that, per class; this is the
% same check made once.
%
% ⚠ `why` names the first disagreement, because "incompatible" with no reason is what sent a
% per-band analysis into a whole-record one in this project before.
%
% See also: rheome.domain.product, rheome.domain.of, rheome.domain.index
%
% Author: Diellor Basha, 2026

    ok = false;  why = "";
    if ~strcmp(dA.kind, dB.kind)
        why = sprintf('kinds differ: %s against %s', dA.kind, dB.kind);  return
    end
    if strcmp(dA.kind, 'product')
        fa = dA.factors;  fb = dB.factors;
        for i = 1:2
            [o, w] = rheome.domain.same(fa{i}, fb{i});
            if ~o, why = sprintf('factor %d: %s', i, w);  return, end
        end
        if dA.nFrames ~= dB.nFrames
            why = sprintf('frame counts differ: %d against %d', dA.nFrames, dB.nFrames);  return
        end
        [o, w] = i_conv(dA.conventions, dB.conventions);
        if ~o, why = w;  return, end
        ok = true;  return
    end
    if ~strcmp(dA.id, dB.id)
        why = sprintf('ids differ: %s against %s', dA.id, dB.id);  return
    end
    ok = true;
end

function [ok, why] = i_conv(ca, cb)
    ok = true;  why = "";
    f = union(fieldnames(ca), fieldnames(cb));
    for i = 1:numel(f)
        ha = isfield(ca, f{i});  hb = isfield(cb, f{i});
        if ha ~= hb
            ok = false;  why = sprintf('convention "%s" is set on one side only', f{i});  return
        end
        if ha && ~isequal(ca.(f{i}), cb.(f{i}))
            ok = false;
            why = sprintf('convention "%s" differs: %s against %s', f{i}, ...
                char(string(ca.(f{i}))), char(string(cb.(f{i}))));
            return
        end
    end
end

% Author: Diellor Basha, 2026
