% DOMAIN  Where a field lives: the domain descriptor and its controlled vocabulary.
%
% A domain says WHAT a carrier is (a triangle mesh, a graph, a point set, an index line),
% HOW BIG it is at each element rank, and WHOSE it is (an id over the topology, a name, the
% file it came from). It carries no coordinates: two fields are on the same domain when the
% ids match, and an operator applies when the kind and the counts do, and neither question
% needs the geometry.
%
%   rheome.domain.kinds     - the four domain kinds and the element ranks each has
%   rheome.domain.of        - a descriptor from a surface / graph / point struct
%   rheome.domain.index     - a descriptor for a line (time, frequency, eigenmodes), with .of
%   rheome.domain.elements  - the element count at a rank
%
% Mirrors the vocabulary of nxr-compute (Manifold, MeshTopology, Domain{vertex,face,edge}),
% so the two sides already agree when a binding exists.
%
% See also: rheome.fieldtype.registry, rheome.operators.registry
%
% Author: Diellor Basha, 2026
