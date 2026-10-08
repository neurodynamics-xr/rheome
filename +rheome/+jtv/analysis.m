function Y = analysis(C, G, info)
% JTV.ANALYSIS  Apply a joint filterbank to a joint spectrum.
%
%   Y = rheome.jtv.analysis(C, G)
%   Y = rheome.jtv.analysis(C, G, info)      % info from rheome.jtv.bank; C may be the struct from rheome.flow.joint
%
% INPUTS:  C     [K x nOmega] joint spectrum, OR the struct from rheome.flow.joint (its .C and .axes are
%                used, so the axes travel with the data rather than beside it)
%          G     [K x nOmega x Nf] filter gains
%          info  (optional) from rheome.jtv.bank. If it carries .axes AND C carries axes, the two are
%                checked with rheome.jtv.compatible before the multiply.
% OUTPUT:  Y     [K x nOmega x Nf] one filtered copy per member
%
% Every member is an elementwise multiply -- the bank costs Nf multiplies and no transforms,
% because both C and G already live on the same (lambda,omega) grid.
%
% ⚠ WHY THE AXES ARE CHECKED AND NOT JUST THE SIZE. Two runs over different bands can retain the
% same number of bins, so [K x nOmega] agreeing proves nothing about them belonging together. So
% can a padded and an unpadded run of the same band -- and those differ in df AND in boundary
% condition. A size check passes all of them silently; rheome.jtv.compatible does not.
%
% See also: rheome.jtv.synthesis, rheome.jtv.dual, rheome.jtv.bank, rheome.jtv.compatible
%
% Author: Diellor Basha, 2026 (after F. Grassi, GSPBox)

    if nargin < 3, info = []; end

    axC = [];
    if isstruct(C)
        if ~isfield(C,'C'), error('jtv:analysis:struct', 'struct input needs a .C field.'); end
        if isfield(C,'axes'), axC = C.axes; end
        C = C.C;
    end

    if ndims(G) < 3, G = reshape(G, size(G,1), size(G,2), []); end
    if ~isequal(size(C), [size(G,1) size(G,2)])
        error('jtv:analysis:size', 'C is %s but the bank is %s.', ...
            mat2str(size(C)), mat2str([size(G,1) size(G,2)]));
    end

    if ~isempty(info) && isstruct(info) && isfield(info,'axes') && ~isempty(info.axes) && ~isempty(axC)
        [ok, why] = rheome.jtv.compatible(info.axes, axC);
        if ~ok
            error('jtv:analysis:axes', ...
                ['the bank was designed on different axes from the spectrum it is applied to:' ...
                 '\n  - %s\nThe sizes match, so nothing else would have caught this.'], ...
                strjoin(why, sprintf('\n  - ')));
        end
    end

    Y = C .* G;
end

% Author: Diellor Basha, 2026
