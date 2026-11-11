classdef Pooling < handle
% CORES.POOLING  Vectors on a folded surface pooled at anchors: carried to the anchor's tangent plane, THEN added (nxr measure.vectors.Pooling).
%
%   P = rheome.cores.Pooling(kit, anchorOf)          % kit: rheome.geom.GeometryKit; anchorOf [N x 1] the anchor vertex of each vertex (0: none)
%   w = P.carry(X, vertices)                         % embedded tangent vectors [k x 3] at vertices -> complex numbers in their anchors' gauges
%   X = P.embed(w, anchors)                          % such numbers (at the anchors) read back as tangent vectors [k x 3] at the anchors (0: zeros)
%
% A tangent vector at a vertex is carried to its anchor (a tile's centre, a place) by the vector heat method's transport
% (rheome.geom.GeometryKit.turnsTo) and only then added, so vectors at different places on a fold are never added in 3-D (a fold turns what is
% parallel along the surface: the raw sum falls short). Defined for any tile size; the transport is along the surface. Vectors of a vertex with
% no anchor, with no angle (zero, or normal to the surface) carry 0.
%
% See also: rheome.geom.GeometryKit, rheome.cores.tilepool, rheome.cores.crossings
%
% Author: Diellor Basha, 2026
    properties (SetAccess = private)
        Kit
        AnchorOf   % [N x 1]
        Turn       % [N x 1] complex: the rotation carrying a vector at each vertex to its anchor
    end
    properties
        Centre = zeros(0, 1)   % [T x 1] a tile pooling's centre vertices (0: empty tile)
        PlaceOf = zeros(0, 1)  % [N x 1] a place pooling's place index of each vertex (0: none)
    end
    methods
        function P = Pooling(kit, anchorOf)
            P.Kit = kit;  P.AnchorOf = anchorOf(:);  P.Turn = kit.turnsTo(P.AnchorOf);
        end
        function w = carry(P, X, vertices)
            vertices = vertices(:);
            w = P.Kit.lift(X, vertices) .* P.Turn(vertices);
            w(~(isfinite(w) & P.AnchorOf(vertices) > 0)) = 0;
        end
        function X = embed(P, w, anchors)
            anchors = anchors(:);  w = w(:);  X = zeros(numel(anchors), 3);  ok = anchors > 0;
            if any(ok), X(ok, :) = P.Kit.embedAtVertices(anchors(ok), w(ok)); end
        end
    end
end
