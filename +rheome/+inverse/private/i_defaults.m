function o = i_defaults(o, Def)
% I_DEFAULTS  Fill unset or empty option fields from a defaults struct.
%
%   o = i_defaults(o, Def)
%
% ⚠ THIS EXISTS AS A FILE BECAUSE rheome.inverse.dirac CALLED IT AND IT WAS NOT IN SCOPE. An identical
% helper is a LOCAL function inside whiten.m, and local functions are file-scoped: whiten.m
% could see it, dirac.m could not. rheome.inverse.dirac therefore failed on EVERY call --
% "Unrecognized function or variable 'i_defaults'" -- with any argument list, because the
% lookup happens before any of its own logic runs. rheome.inverse.mne is unaffected: it fills its
% defaults with an inline loop instead.
%
% A function in +inverse/private/ is visible to every function in +inverse/, which is the scope
% this needs. whiten.m keeps its local copy, which shadows this one inside that file.
%
% See also: rheome.inverse.dirac, rheome.inverse.mne
%
% Author: Diellor Basha, 2026

    f = fieldnames(Def);
    for i = 1:numel(f)
        if ~isfield(o, f{i}) || isempty(o.(f{i}))
            o.(f{i}) = Def.(f{i});
        end
    end
end
% Author: Diellor Basha, 2026
