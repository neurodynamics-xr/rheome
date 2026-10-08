classdef tFiltersSchrodinger < matlab.unittest.TestCase
    % rheome.filters.schrodinger: the unitary member of the propagator family.
    %
    % Author: Diellor Basha, 2026

    properties (Constant)
        lam = linspace(0, 4, 40)'
        t   = linspace(0, 3, 60)
    end

    methods (Test)

        function itIsUnitary(tc)
            % ⭐ the property that distinguishes it: it turns, it does not spread
            G = rheome.filters.schrodinger(tc.lam, tc.t);
            tc.verifySize(G, [numel(tc.lam) numel(tc.t)]);
            tc.verifyEqual(abs(G), ones(size(G)), 'AbsTol', 1e-12);
            tc.verifyFalse(isreal(G));
        end

        function diffusionIsTheSameOperatorInRealTime(tc)
            % exp(-lambda t) against exp(-i lambda t): the Wick rotation, and nothing else
            H = rheome.filters.diffusion(tc.lam, tc.t);
            G = rheome.filters.schrodinger(tc.lam, tc.t);
            tc.verifyEqual(abs(log(max(H,realmin))), abs(angle(G)), 'AbsTol', 1e-9);
        end

        function eachModeTurnsAtItsOwnEigenvalue(tc)
            % ⭐ the dispersion relation omega = lambda, which is what makes it a straight
            %   ridge in the joint plane
            G = rheome.filters.schrodinger(tc.lam, tc.t, 1, 1);      % lmax = 1, so omega = lambda
            for k = [5 17 33]
                ph = unwrap(angle(G(k,:)));
                rate = -median(diff(ph))/median(diff(tc.t));
                tc.verifyEqual(rate, tc.lam(k), 'RelTol', 1e-6, 'AbsTol', 1e-9);
            end
        end

        function zeroLagIsTheIdentity(tc)
            G = rheome.filters.schrodinger(tc.lam, 0);
            tc.verifyEqual(G, ones(numel(tc.lam),1), 'AbsTol', 1e-12);
        end

        function hbarSetsTheRate(tc)
            G1 = rheome.filters.schrodinger(tc.lam, tc.t, 1);
            G2 = rheome.filters.schrodinger(tc.lam, tc.t, 2);
            tc.verifyEqual(angle(G2(10,:)), angle(G1(10,:))/2, 'AbsTol', 1e-9);
            tc.verifyError(@() rheome.filters.schrodinger(tc.lam, tc.t, 0), ?MException);
        end

    end
end

% Author: Diellor Basha, 2026
