function [ok, why] = check(id, dom, varargin)
% SELECTION.CHECK  Is this selection applicable to this domain, and may a coarser cell be summed?
%
%   [ok, why] = rheome.selection.check("cortex_tile", e.domains.source)
%   [ok, why] = rheome.selection.check("graph_wavelet", d, Merging=true)
%
% ⭐ THE SECOND QUESTION IS THE USEFUL ONE. Asking whether a selection FITS a domain catches a
% mismatched axis; asking whether it MERGES catches the mistake that actually happens -- summing a
% coarser cell from finer ones under a kernel, where the frame bound is not 1 and the sum is not the
% whole. Pass Merging=true to require it.
%
% ⚠ A true answer here permits merging; it does not make every measurement mergeable. Only
% quantities linear in the accumulated statistic roll up even under a partition (rheome.flow.windowtable:
% `energy` does, peak and crest and vortex count do not).
%
% Author: Diellor Basha, 2026
    p = inputParser;
    p.addParameter('Merging', false, @islogical);
    p.parse(varargin{:});
    s = rheome.selection.describe(id);
    ok = true;  why = "";
    % ⭐ a product selection cuts a product domain; a single-axis one may also be applied to a
    % product's matching factor, which is what makes a spatial tile usable on (cortex x time).
    if isstruct(dom) && isfield(dom, 'kind') && strcmp(dom.kind, 'product') ...
            && ~strcmp(s.domain_kind, 'product')
        dom = dom.factors{1};
    end
    if isstruct(dom) && isfield(dom, 'kind') && ~strcmp(dom.kind, s.domain_kind)
        ok = false;
        why = sprintf('%s cuts a "%s" domain, not "%s"', id, s.domain_kind, dom.kind);
        return
    end
    if p.Results.Merging && ~s.merges
        ok = false;
        why = sprintf('%s is a %s: a coarser cell is NOT the sum of finer ones', id, s.kind);
    end
end

% Author: Diellor Basha, 2026
