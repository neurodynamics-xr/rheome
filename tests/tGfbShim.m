classdef tGfbShim < matlab.unittest.TestCase
% rheome.filters.frame now delegates to graphfilterbank. These tests prove the conversion
% changed NOTHING, field for field, against tests/gfbGoldenFrame.mat -- captured from
% the pre-shim implementation.

    properties
        Golden
    end

    methods (TestClassSetup)
        function loadGolden(tc)
            S = load(fullfile(fileparts(which('gfbFixture')), 'gfbGoldenFrame.mat'));
            tc.Golden = S.G;
        end
    end

    methods (Test)

        function everyGoldenCaseMatchesFieldForField(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            G  = tc.Golden;
            numFields = {'Lrange','Centers','t','Sigma','Gamma','SigmaTrue', ...
                         'MassLost','SigmaFloor','SigmaMax','LminEff','PerOctave'};
            for i = 1:numel(G)
                g = G(i);
                lbl = sprintf('case %d (%s)', i, g.family);
                f = rheome.filters.frame(g.family, g.Nf, g.lrange, g.opts{:}, 'warn', false);

                tc.verifyEqual(numel(f.g), g.M,       ['member count, ' lbl]);
                tc.verifyEqual(f.Nf,       g.NfField, ['Nf field, ' lbl]);
                tc.verifyEqual(f.Family,   g.family,  ['Family, ' lbl]);
                tc.verifyEqual(f.Usable,   g.Usable,  ['Usable, ' lbl]);

                % The gains are the load-bearing check: the handles must evaluate
                % identically, bit for bit.
                tc.verifyEqual(rheome.filters.frame_gains(f, lam), g.H, 'AbsTol', 0, ...
                    ['gains, ' lbl]);

                for fld = numFields
                    tc.verifyEqual(f.(fld{1}), g.(fld{1}), 'RelTol', 1e-12, ...
                        sprintf('%s, %s', fld{1}, lbl));
                end
            end
        end

        function truncatedWarningStillFires(tc)
            % sphere_validate.m asserts on this exact identifier.
            fx = gfbFixture();
            tc.verifyWarning(@() rheome.filters.frame('mexhat', 7, fx.Lambda), ...
                'filters:frame:truncated');
        end

        function sparseWarningStillFires(tc)
            % omega_rerun.m and flow_curl_methods.m suppress this exact identifier.
            fx = gfbFixture();
            tc.verifyWarning(@() rheome.filters.frame('mexhat', 3, fx.Lambda), ...
                'filters:frame:sparse');
        end

        function warnFalseIsSilent(tc)
            fx = gfbFixture();
            tc.verifyWarningFree(@() rheome.filters.frame('mexhat', 3, fx.Lambda, 'warn', false));
        end

        function errorIdentifiersArePreserved(tc)
            fx = gfbFixture();
            tc.verifyError(@() rheome.filters.frame('nonsense', 5, fx.Lambda), 'filters:frame:family');
            tc.verifyError(@() rheome.filters.frame('mexhat', 5, fx.Lambda, 'tmin','bogus'), 'filters:frame:tmin');
            tc.verifyError(@() rheome.filters.frame('mexhat', 5, []), 'filters:frame:lrange');
            tc.verifyError(@() rheome.filters.frame('mexhat', 5, [100 2]), 'filters:frame:lrange');
        end

        function downstreamConsumersStillWork(tc)
            % frame_gains/bounds/analysis/synthesis/scalogram consume the STRUCT, so they
            % should need no change. Assert that end to end.
            fx = gfbFixture();  lam = fx.Lambda;
            b  = struct('Phi',fx.Phi, 'Lambda',lam, 'Mass',fx.Mass, 'nV',fx.nV);
            f  = rheome.filters.frame('mexhat', 6, lam, 'warn', false);
            % Band-limited: the basis retains only K modes, so a full-rank random field
            % has content the transform cannot represent and 'dual' correctly returns its
            % projection rather than F.
            F  = fx.Phi * randn(numel(lam), 4);
            W  = rheome.filters.frame_analysis(b, F, f);
            tc.verifySize(W, [fx.nV, 4, numel(f.g)]);
            tc.verifyEqual(rheome.filters.frame_synthesis(b, W, f, 'dual'), F, 'AbsTol', 1e-8);
            tc.verifyNotEmpty(rheome.filters.frame_bounds(f, lam).A);
            tc.verifyNotEmpty(rheome.filters.frame_scalogram(b, F, f).energy);
            tc.verifySize(rheome.filters.localize(b, [1;2], rheome.filters.frame_gains(f, lam)), ...
                [fx.nV, 2, numel(f.g)]);
        end

    end
end

% Author: Diellor Basha, 2026
