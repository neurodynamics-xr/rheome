classdef tCrosstalk < matlab.unittest.TestCase
% rheome.flow.crosstalk -- where an ROI's measured flow actually comes from.
% Author: Diellor Basha, 2026

    methods (Test)

        function aPerfectlyResolvedRegionIsFullyOwn(tc)
            % If the kernel reads each vertex from its own source only, every ROI is 100% own.
            V=12; C=V; P=logical([ones(1,6) zeros(1,6); zeros(1,6) ones(1,6)]);
            K=eye(V); Gain=zeros(C,3*V);
            for v=1:V, Gain(v,3*(v-1)+1)=1; end
            [own,dom]=rheome.flow.crosstalk(K,Gain,P);
            tc.verifyEqual(own, [1;1], 'AbsTol', 1e-12);
            tc.verifyEqual(dom, [1;2]);
        end

        function aRegionReadingAnotherRegionsSourcesIsFlaggedAsSuch(tc)
            % The failure this exists to catch: ROI 1's measurement is entirely driven by
            % sources in ROI 2. own -> 0, and `dominant` names the true origin.
            V=12; P=logical([ones(1,6) zeros(1,6); zeros(1,6) ones(1,6)]);
            K=eye(V); Gain=zeros(V,3*V);
            for v=1:6,  Gain(v,3*(v+5)+1)=1; end     % ROI 1 rows read ROI 2 sources
            for v=7:12, Gain(v,3*(v-1)+1)=1; end
            [own,dom]=rheome.flow.crosstalk(K,Gain,P);
            tc.verifyLessThan(own(1), 1e-9);
            tc.verifyEqual(dom(1), 2);
            tc.verifyEqual(own(2), 1, 'AbsTol', 1e-12);
        end

        function ownFractionIsBoundedInZeroOne(tc)
            rng(2); V=30; C=8;
            P=false(3,V); P(1,1:10)=true; P(2,11:20)=true; P(3,21:30)=true;
            own=rheome.flow.crosstalk(randn(V,C), randn(C,3*V), P);
            tc.verifyTrue(all(own>=0 & own<=1));
        end

        function anEmptyRegionGivesNaNRatherThanZero(tc)
            % Zero would read as "fully leakage"; the truth is "no measurement".
            V=10; P=false(2,V); P(1,1:5)=true;
            own=rheome.flow.crosstalk(randn(V,4), randn(4,3*V), P);
            tc.verifyTrue(isnan(own(2)));
            tc.verifyFalse(isnan(own(1)));
        end

        function rejectsAMismatchedLeadfield(tc)
            tc.verifyError(@() rheome.flow.crosstalk(randn(10,4), randn(5,30), true(1,10)), ...
                'flow:crosstalk:size');
        end

    end
end

% Author: Diellor Basha, 2026
