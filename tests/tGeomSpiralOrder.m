classdef tGeomSpiralOrder < matlab.unittest.TestCase
% rheome.geom.spiralorder: points on the spiral come back at their own parameter, the order runs pole to
% pole, it follows the pole axis, and it sorts tiles as it sorts vertices.
%
% Author: Diellor Basha, 2026
    methods (Test)
        function pointsOnTheSpiralComeBackAtTheirParameter(t)
            T = 12;  s = linspace(0.02, 0.98, 500)';
            for k = ["colat" "area" "loxodrome"]
                switch k
                    case "colat", th = pi*s;
                    case "area",  th = acos(1 - 2*s);
                    case "loxodrome"
                        y0 = asinh(cotd(2));  th = acot(sinh(y0*(1 - 2*s)));  th(th < 0) = th(th < 0) + pi;
                end
                ph = 2*pi*T*s;
                X = [sin(th).*cos(ph) sin(th).*sin(ph) cos(th)];
                [o, tt] = rheome.geom.spiralorder(3*X, Turns=T, Kind=k);
                t.verifyEqual(tt, s, 'AbsTol', 1e-9, char(k));
                t.verifyEqual(o, (1:500)', char(k));
            end
        end
        function itRunsFromNorthToSouthOnAMesh(t)
            [V, ~] = rheome.geom.icosphere(4);
            [o, tt, w] = rheome.geom.spiralorder(V, Turns=8);
            t.verifyEqual(sort(o), (1:size(V,1))');
            t.verifyTrue(issorted(tt(o)));
            t.verifyGreaterThan(V(o(1),3), 0.99);  t.verifyLessThan(V(o(end),3), -0.99);
            t.verifyGreaterThanOrEqual(min(w), -1);  t.verifyLessThanOrEqual(max(w), 8);
            % a winding is a strip: its colatitudes stay within one gap (180/8 deg) of its centre
            c = acosd(V(:,3));  in = w >= 1 & w <= 6;
            t.verifyLessThan(max(abs(c(in) - 180*(w(in) + mod(atan2(V(in,2),V(in,1)),2*pi)/(2*pi))/8)), 180/8/2 + 1e-9);
        end
        function itFollowsThePoleAxis(t)
            [V, ~] = rheome.geom.icosphere(3);
            R = [0 0 1; 0 -1 0; 1 0 0];                     % proper rotation swapping z and x: it carries
            %                                                 the z chart's longitude origin (x) to the x chart's (z)
            [oz, tz] = rheome.geom.spiralorder(V, Turns=6);
            [ox, tx] = rheome.geom.spiralorder(V*R', Turns=6, Axis=[1 0 0]);
            t.verifyEqual(tx, tz, 'AbsTol', 1e-12);
            t.verifyEqual(ox, oz);
        end
        function tilesSortLikeTheirVertices(t)
            [V, ~] = rheome.geom.icosphere(5);
            P = rheome.geom.spherepatches(V, Level=2);
            o = rheome.geom.spiralorder(P.centre, Turns=4);
            [~, tv] = rheome.geom.spiralorder(V, Turns=4);
            tm = accumarray(P.patch, tv, [], @median);      % a tile's median vertex position
            % 0.98, not 1: a tile straddling a winding holds vertices one winding apart (measured 0.990)
            t.verifyGreaterThan(corr(tm(o), (1:numel(o))', 'Type', 'Spearman'), 0.98);
        end
    end
end

% Author: Diellor Basha, 2026
