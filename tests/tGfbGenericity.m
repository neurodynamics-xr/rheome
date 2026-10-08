classdef tGfbGenericity < matlab.unittest.TestCase
% The point of the exercise: the class must work with none of this repo present.

    methods (Test)

        function classReferencesNoRepoPackage(tc)
            % Mechanical check. @graphfilterbank and +graphtransform must be portable.
            root = fileparts(fileparts(which('gfbFixture')));
            pkgs = {'filters\.', 'eigen\.', 'flow\.', 'operators\.', 'detect\.', ...
                    'forward\.', 'inverse\.', 'source\.', 'jtv\.', 'geom\.', ...
                    'load\.', 'import\.', 'show\.', 'spectral\.', 'connectome\.'};
            files = [dir(fullfile(root, '+rheome', '@graphfilterbank', '**', '*.m')); ...
                     dir(fullfile(root, '+rheome', '@jointfilterbank', '**', '*.m')); ...
                     dir(fullfile(root, '+rheome', '+graphtransform', '*.m'))];
            tc.verifyNotEmpty(files);
            for f = files'
                src = fileread(fullfile(f.folder, f.name));
                src = regexprep(src, '^\s*%.*$', '', 'lineanchors');   % strip comments
                for pk = pkgs
                    tc.verifyEmpty(regexp(src, pk{1}, 'once'), ...
                        sprintf('%s references %s', f.name, pk{1}));
                end
            end
        end

        function combinatorialGraphLaplacianWorks(tc)
            % No mesh, no mass matrix, dimensionless lambda.
            rng(0);
            n = 60;
            A = double(rand(n) < 0.08);  A = triu(A,1);  A = A + A';
            L = diag(sum(A,2)) - A;
            [U, D] = eig(L);
            lam = diag(D);
            T   = rheome.graphtransform.eigen(U, speye(n), lam);
            g   = rheome.graphfilterbank(lam, 'Wavelet','itersine', 'NumFilters',8, 'Transform',T);
            x   = randn(n, 1);
            tc.verifyEqual(iwt(g, wt(g, x), 'dual'), x, 'AbsTol', 1e-8);
        end

        function designMethodsWorkWithNoTransform(tc)
            g = rheome.graphfilterbank(100, 'Wavelet','mexhat');
            tc.verifyNotEmpty(graphfilters(g));
            tc.verifyNotEmpty(framebounds(g));
            tc.verifyNotEmpty(scales(g));
            tc.verifyNotEmpty(widths(g));
            tc.verifyNotEmpty(centerWavenumbers(g));
            tc.verifyNotEmpty(wavelengths(g));
            tc.verifyNotEmpty(qfactor(g));
            tc.verifyNotEmpty(powerbw(g));
            tc.verifyNotEmpty(spectralResponse(g));
            tc.verifyWarningFree(@() isframetight(g));
        end

        function quaternionFieldFiltersThroughNormOnly(tc)
            % Dirac differs from LBO ONLY in Transform.norm.
            rng(0);
            nV = 30;  K = 12;
            Phi = orth(randn(4*nV, K));
            lam = linspace(1, 80, K)';
            T   = rheome.graphtransform.eigen(Phi, speye(4*nV), lam, 'nV', nV);
            g   = rheome.graphfilterbank(lam, 'Wavelet','itersine', 'NumFilters',6, 'Transform',T);
            F   = Phi * randn(K, 3);
            tc.verifyEqual(iwt(g, wt(g, F), 'dual'), F, 'AbsTol', 1e-8);
            tc.verifySize(vertexSpectrum(g, F), [g.NumMembers, 3]);
            tc.verifySize(scaleSpectrum(g, F),  [nV, 3]);
        end

        function complexFieldFiltersUnchanged(tc)
            % The connection Laplacian case: complex Phi, no other difference.
            rng(0);
            n = 40;  K = 10;
            Phi = orth(randn(n,K) + 1i*randn(n,K));
            lam = linspace(1, 50, K)';
            T   = rheome.graphtransform.eigen(Phi, speye(n), lam);
            g   = rheome.graphfilterbank(lam, 'Wavelet','itersine', 'NumFilters',6, 'Transform',T);
            F   = Phi * (randn(K,2) + 1i*randn(K,2));
            tc.verifyEqual(iwt(g, wt(g, F), 'dual'), F, 'AbsTol', 1e-8);
        end

        function jointBankWorksOnAPlainGraphLaplacian(tc)
            % No mesh, no mass matrix, no repo package anywhere.
            rng(0);
            n = 50;
            A = double(rand(n) < 0.1);  A = triu(A,1);  A = A + A';
            L = diag(sum(A,2)) - A;
            [U, D] = eig(L);  lam = diag(D);
            gfb = rheome.graphfilterbank(lam, 'Wavelet','itersine', 'NumFilters',4, ...
                                  'Transform', rheome.graphtransform.eigen(U, speye(n), lam));
            j = rheome.jointfilterbank(gfb, 'SignalLength', 32, 'SamplingFrequency', 32);
            C = randn(n, j.NumOmega);
            tc.verifyEqual(iwt(j, wt(j, C), 'dual'), C, 'AbsTol', 1e-8);
        end

        function filtersFrameShimStillWorks(tc)
            % Every existing caller must keep working, unchanged.
            fx = gfbFixture();  lam = fx.Lambda;
            f  = rheome.filters.frame('mexhat', 7, lam, 'warn', false);
            tc.verifyEqual(numel(f.g), 8);
            tc.verifySize(rheome.filters.frame_gains(f, lam), [numel(lam), 8]);
            b = rheome.filters.frame_bounds(f, lam);
            tc.verifyTrue(isfield(b,'A') && isfield(b,'B') && isfield(b,'Tightness'));
            tc.verifyTrue(isfield(f,'Sigma') && isfield(f,'MassLost') && isfield(f,'Usable'));
        end

    end
end

% Author: Diellor Basha, 2026
