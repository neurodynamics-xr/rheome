classdef tJfbCompose < matlab.unittest.TestCase

    methods (Test)

        function memberCountIsTheProductOfTheThreeLists(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Bands', [4 8; 8 13], ...
                                 'JointKernels', {@(l,w) 1, @(l,w) 1});
            tc.verifyEqual(j.NumGraph, fx.gfb.NumMembers);
            tc.verifyEqual(j.NumTime,  2);
            tc.verifyEqual(j.NumJoint, 2);
            tc.verifyEqual(j.NumMembers, fx.gfb.NumMembers * 2 * 2);
        end

        function indexMapsEveryMemberToItsTriple(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[4 8; 8 13]);
            tc.verifySize(j.Index, [j.NumMembers, 3]);
            tc.verifyEqual(sort(unique(j.Index(:,1))).', 1:j.NumGraph);
            tc.verifyEqual(sort(unique(j.Index(:,2))).', 1:j.NumTime);
            tc.verifyEqual(size(unique(j.Index, 'rows'), 1), j.NumMembers);
        end

        function labelsAreUniqueAndCountMatches(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[4 8; 8 13]);
            tc.verifyNumElements(j.Labels, j.NumMembers);
            tc.verifyNumElements(unique(j.Labels), j.NumMembers);
        end

        function defaultsToOneAllPassTimeMemberAndOneUnitKernel(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength', 64);
            tc.verifyEqual(j.NumTime,  1);
            tc.verifyEqual(j.NumJoint, 1);
            tc.verifyEqual(j.NumMembers, fx.gfb.NumMembers);
            tc.verifyTrue(j.Separable);
        end

        function bandsBecomeRaisedCosineTimeMembers(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Bands', [8 13]);
            tc.verifyEqual(j.NumTime, 1);
            % TimeFilters is the raw input; TimeMembers is the RESOLVED list.
            h = j.TimeMembers{1};
            tc.verifyEqual(h(2*pi*8),    0, 'AbsTol', 0);
            tc.verifyEqual(h(2*pi*10.5), 1, 'RelTol', 1e-12);
        end

        function explicitTimeFiltersCoexistWithBands(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                     'Bands', [8 13], 'TimeFilters', {@(w) ones(size(w))});
            tc.verifyEqual(j.NumTime, 2);
        end

        function separableIsFalseWithANonUnitKernel(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, ...
                     'JointKernels', {@(l,w) exp(-(w - 0.8*sqrt(l)).^2 / (2*0.2^2))});
            tc.verifyFalse(j.Separable);
        end

        function aSingleHandleIsAcceptedWithoutACell(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, ...
                     'TimeFilters', @(w) ones(size(w)), 'JointKernels', @(l,w) 1);
            tc.verifyEqual(j.NumTime, 1);
            tc.verifyEqual(j.NumJoint, 1);
        end

        function labelCountMismatchIsRefused(tc)
            fx = jfbFixture();
            tc.verifyError(@() rheome.jointfilterbank(fx.gfb, 'SignalLength',64, ...
                'Labels', {'only-one'}), 'jointfilterbank:labels');
        end

    end
end

% Author: Diellor Basha, 2026
