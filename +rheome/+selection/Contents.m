% SELECTION  How we cut a manifold: tiles, windows, bands, patches, ROIs and wavelet members as
% ONE type -- a weight over the elements of a domain.
%
% The fourth pillar beside +domain (where a field lives), +fieldtype (what it is) and
% rheome.operators.registry (which arrows move it).
%
%   rheome.selection.registry  - every selection, with the property that matters: does it MERGE
%   rheome.selection.describe  - one selection by id
%   rheome.selection.check     - does it fit this domain, and may a coarser cell be summed from finer ones
%
% ⭐ partition merges exactly (rheome.geom.tree, dyadic tiles, constant-Q bands, Weyl octaves);
%   kernel does not (a frame you invert, not a sum -- one joint wavelet member keeps 1/21 of its
%   cell); subset does not (an ROI list is neither disjoint nor exhaustive).
%
% See also: rheome.domain.kinds, rheome.fieldtype.registry, rheome.operators.registry
%
% Author: Diellor Basha, 2026
