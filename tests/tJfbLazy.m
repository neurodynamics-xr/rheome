classdef tJfbLazy < matlab.unittest.TestCase
% The central property: a member is evaluated on demand from its factors.

    methods (Test)

        function separableMemberIsExactlyTheOuterProduct(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                                 'Bands', [8 13]);
            m  = 2;
            ig = j.Index(m,1);  it = j.Index(m,2);
            gv = gain(fx.gfb, ig);
            hv = j.TimeMembers{it};
            expected = gv(j.Lambda) * hv(j.Omega);
            tc.verifyEqual(jointfilters(j, m), expected, 'AbsTol', 0);
        end

        function memberShapeIsKByNumOmega(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[4 8; 8 13]);
            tc.verifySize(jointfilters(j, 1), [numel(j.Lambda), j.NumOmega]);
        end

        function nonSeparableKernelChangesTheMember(tc)
            fx = jfbFixture();
            sep = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64);
            spd = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'SamplingFrequency',64, ...
                     'JointKernels', {@(l,w) exp(-(w - 5*sqrt(l)).^2 / (2*10^2))});
            tc.verifyNotEqual(jointfilters(spd,1), jointfilters(sep,1));
        end

        function aScalarKernelBroadcasts(tc)
            % @(l,w) 1 returns a scalar; it must expand, not error.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'JointKernels', {@(l,w) 1});
            tc.verifySize(jointfilters(j, 1), [numel(j.Lambda), j.NumOmega]);
        end

        function aKernelOfTheWrongSizeIsRefused(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, ...
                                 'JointKernels', {@(l,w) ones(3,4)});
            tc.verifyError(@() jointfilters(j, 1), 'jointfilterbank:kernelSize');
        end

        function memberIndexIsValidated(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64);
            tc.verifyError(@() jointfilters(j, 0),               'jointfilterbank:member');
            tc.verifyError(@() jointfilters(j, j.NumMembers+1),  'jointfilterbank:member');
        end

        function fullMaterialisationWorksWhenSmall(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[8 13]);
            G  = jointfilters(j);
            tc.verifySize(G, [numel(j.Lambda), j.NumOmega, j.NumMembers]);
            for m = 1:j.NumMembers
                tc.verifyEqual(G(:,:,m), jointfilters(j,m), 'AbsTol', 0);
            end
        end

        function theBankIsFullyUsableBelowTheGuard(tc)
            % THE central property of the design: with MaxBytes set so small that the
            % product cannot be formed, everything still works EXCEPT forming it.
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'Bands',[8 13], ...
                                 'MaxBytes', 1000);
            tc.verifyNotEmpty(jointfilters(j, 1));          % one member: fine
            tc.verifyNotEmpty(framebounds(j));              % reductions: fine
            tc.verifyError(@() jointfilters(j), 'jointfilterbank:tooLarge');
        end

        function theGuardNamesTheSize(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'MaxBytes', 1000);
            try
                jointfilters(j);
                tc.verifyFail('expected jointfilterbank:tooLarge');
            catch ME
                tc.verifyEqual(ME.identifier, 'jointfilterbank:tooLarge');
                tc.verifySubstring(ME.message, 'MaxBytes');
                tc.verifySubstring(ME.message, 'Force');
            end
        end

        function forceOverridesTheGuard(tc)
            fx = jfbFixture();
            j  = rheome.jointfilterbank(fx.gfb, 'SignalLength',64, 'MaxBytes', 1000);
            tc.verifyNotEmpty(jointfilters(j, [], 'Force', true));
        end

    end
end

% Author: Diellor Basha, 2026
