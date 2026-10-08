function d = product(dA, dB, varargin)
% DOMAIN.PRODUCT  A joint manifold: a domain crossed with an index line.
%
%   d = rheome.domain.product(dCortex, rheome.domain.index(nT, Name="time"))
%   d = rheome.domain.product(dModes, dFreq, Conventions=struct('half','positive','boundary','ring'))
%
% ⭐ WHY THIS EXISTS. A joint plane -- (cortex x time), (eigenmode x frequency), (time x frequency) --
% was represented four different and unconnected ways: by juxtaposition in +pipeline (a spatial domain
% plus a Frames axis carried alongside), by a conjugacy pointer in rheome.domain.index(Of=), by a class-local
% axes struct in @jointfilterbank with its own equality check, and by a composite key in the tile
% store. None of them was a TYPE, so "are these two analyses on the same plane" had a different answer
% in each place. This makes it one descriptor and one equality (rheome.domain.same).
%
% ⭐ THE DESIGN CONSTRAINT THAT MADE IT CHEAP. A product mirrors FACTOR 1's element counts in
% .nV/.nE/.nF, so every existing shape lambda in rheome.operators.registry -- @(d)[d.nV,d.nV] and the other
% 11 -- keeps working unchanged on a product domain. The product count is reached through
% rheome.domain.elements(d, 'element') and .nFrames, which nothing older asks for. Adding the kind therefore
% breaks no arrow: 20 of the 34 are within one manifold and act on the spatial factor, identity on the
% other, and they cannot tell the difference.
%
% ⚠ FACTOR 2 MUST BE AN INDEX LINE. Every joint plane in this project is (something x a line): time,
% frequency or modes. A product of two complexes is not supported and errors rather than half-working.
%
% ⚠⚠ CONVENTIONS ARE PART OF THE DOMAIN, and that is the point rather than an extra. Two banks can
% agree on every count and still be incompatible because one holds the positive half of the spectrum
% and the other the whole of it, or one treats time as a ring and the other as a path.
% @jointfilterbank/iscompatible checks .half and .boundary "because a mismatch there is silent and
% fatal" -- that check now lives in the type, so it is made once and not per class.
%
% INPUTS
%   dA, dB       domain descriptors; dB.kind must be 'index'
%   Conventions  (optional) struct of convention fields compared by rheome.domain.same
%   Name         (optional) a label
%
% OUTPUT
%   d  .kind='product'  .factors={dA,dB}  .nV/.nE/.nF mirroring dA  .nFrames=dB.nV
%      .conventions  .name  .id (hashed from both factors and the conventions)
%
% See also: rheome.domain.same, rheome.domain.elements, rheome.domain.index, rheome.domain.kinds, rheome.selection.registry
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Conventions', struct(), @isstruct);
    p.addParameter('Name', "", @(x) isstring(x) || ischar(x));
    p.parse(varargin{:});
    o = p.Results;

    if ~isstruct(dA) || ~isfield(dA,'kind'), error('domain:product:factor','dA is not a domain.'); end
    if ~isstruct(dB) || ~isfield(dB,'kind'), error('domain:product:factor','dB is not a domain.'); end
    if ~strcmp(dB.kind, 'index')
        error('domain:product:second', ...
            'the second factor must be an index line (time, frequency, modes), not a %s.', dB.kind);
    end

    nm = char(o.Name);
    if isempty(nm)
        a = dA.name;  if isempty(a), a = dA.kind; end
        b = dB.name;  if isempty(b), b = 'index'; end
        nm = sprintf('%s x %s', a, b);
    end
    d = struct('kind', 'product', 'factors', {{dA, dB}}, ...
               'nV', dA.nV, 'nE', dA.nE, 'nF', dA.nF, 'nFrames', dB.nV, ...
               'id', '', 'name', nm, 'source', '', 'units', '', 'of', '', ...
               'conventions', o.Conventions);
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(sprintf('product|%s|%s|%s', dA.id, dB.id, i_convstr(o.Conventions))));
    dg = typecast(md.digest(), 'uint8');
    d.id = lower(reshape(dec2hex(dg(1:4), 2)', 1, []));
end

function s = i_convstr(c)
    f = sort(fieldnames(c));  parts = cell(1, numel(f));
    for i = 1:numel(f)
        v = c.(f{i});
        if isnumeric(v) || islogical(v), v = mat2str(v); else, v = char(string(v)); end
        parts{i} = sprintf('%s=%s', f{i}, v);
    end
    s = strjoin(parts, ';');
end

% Author: Diellor Basha, 2026
