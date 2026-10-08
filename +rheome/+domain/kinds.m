function T = kinds()
% DOMAIN.KINDS  The controlled vocabulary of domains, and which element ranks each has.
%
%   T = rheome.domain.kinds()
%
% A domain is WHERE a field lives. This repo has four, and they differ in which element
% ranks exist on them -- which is the whole reason an operator's arrow can be checked:
%
%   complex  a triangle mesh: vertices, edges and faces (the cortex, a sphere)
%   graph    a weighted graph: vertices and edges, no faces (the connectome)
%   points   an unstructured point set: vertices only (a sensor array, before a graph)
%   index    a one-dimensional index set: elements only (time samples, frequency bins,
%            eigenmode indices -- the "line" a transform lands on)
%   product  a joint manifold: one of the above crossed with an index line (cortex x time,
%            eigenmode x frequency). ⭐ It MIRRORS FACTOR 1's element counts, so every existing
%            shape lambda works on it unchanged; the product count is rheome.domain.elements(d,'element')
%            and .nFrames. ⚠ factor 2 must be an index line -- see rheome.domain.product.
%
% ⚠ A PRODUCT'S RANKS ARE ITS FIRST FACTOR'S, PLUS 'element' for the product itself, so the ranks
% listed here are the union and rheome.domain.elements defers the real check to factor 1. A
% (index x index) product has no 'vertex' rank even though the table lists one.
%
% ⚠ A SENSOR ARRAY IS USUALLY A COMPLEX HERE. rheome.sensors.graph returns .Vertices/.Faces/.nV
% exactly as the mesh code expects, which is why the flow operators accept it with no
% adapter; so a sensor array is 'points' only before that triangulation is built.
%
% See also: rheome.domain.of, rheome.domain.index, rheome.fieldtype.registry
%
% Author: Diellor Basha, 2026

    r = { 'complex', {'vertex','edge','face'}, 'triangle mesh (simplicial 2-complex)'
          'graph',   {'vertex','edge'},        'weighted graph, no faces'
          'points',  {'vertex'},               'unstructured point set'
          'index',   {'element'},              'one-dimensional index set (time, frequency, modes)'
          'product', {'vertex','edge','face','element'}, 'a domain crossed with an index line (see rheome.domain.product)' };
    T = table(string(r(:,1)), r(:,2), string(r(:,3)), ...
              'VariableNames', {'kind','ranks','description'});
end
% Author: Diellor Basha, 2026
