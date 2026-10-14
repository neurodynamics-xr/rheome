classdef tAtlasTime < matlab.unittest.TestCase
% The atlas's time tiles and its alpha wavelet, as detection uses them.
    methods (Test)
        function alphaTile(tc)
            A = rheome.detect.atlastime(10);
            tc.verifyEqual(A.level, -19);  tc.verifyEqual(A.P, 0.164794921875, AbsTol=1e-15);
            tc.verifyEqual(rheome.detect.atlastime(sqrt(8*16)).level, -19);       % MS1's band, same tile
            tc.verifyEqual(rheome.detect.atlastime(A.band(2)).level, -20);         % edges are half-open
        end
        function tilesAverageAndAlign(tc)
            P = 0.1648;  r = 300;  y = 1:1000;
            [Yt, i] = rheome.detect.timetiles(y, r, P, T0=0.5);
            tc.verifyEqual(i.k(1), ceil(0.5/P - 0.5/(r*P)));
            t = 0.5 + (0:999)/r;  k1 = floor(t/P + 1e-9) == i.k(1);
            tc.verifyEqual(Yt(1), mean(y(k1)), AbsTol=1e-12);
            tc.verifyLessThanOrEqual(abs(diff(i.mid(1:2)) - P*r), 1);
        end
        function alphaWavelet(tc)
            fs = 600;  t = (0:6000-1)/fs;  A = rheome.detect.atlastime(10);
            xa = rheome.filters.atlasband(3*cos(2*pi*A.fc*t), fs, 10);
            tc.verifyEqual(median(abs(xa)), 3, RelTol=1e-3);                       % amplitude kept at the centre
            g = @(f) median(abs(rheome.filters.atlasband(cos(2*pi*f*t), fs, 10)));
            tc.verifyEqual([g(A.band(1)) g(A.band(2))], [1 1]/sqrt(2), AbsTol=0.02);  % half power at 1/P, 2/P
            tc.verifyLessThan(g(30), 1e-6);
        end
    end
end
