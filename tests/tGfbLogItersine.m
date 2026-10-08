classdef tGfbLogItersine < matlab.unittest.TestCase
% The log-itersine family of graphfilterbank: tight on the whole spectrum, members
% spaced by voices in log wavenumber, edge members completing the frame -- the temporal
% bank's construction on a graph spectrum.
%
% Author: Diellor Basha, 2026

    methods (Test)

        function itIsTightOnASuppliedSpectrum(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            for V = [1 2 3]
                g = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', V);
                b = framebounds(g, lam);
                tc.verifyEqual(b.A, 1, 'AbsTol', 1e-10, sprintf('V=%d', V));
                tc.verifyEqual(b.B, 1, 'AbsTol', 1e-10, sprintf('V=%d', V));
                tc.verifyTrue(isframetight(g));
                tc.verifyEqual(g.MassLost, zeros(1, g.NumMembers));
            end
        end

        function membersAreLogSpacedAndVoicesAreHonoured(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g1 = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', 1);
            g2 = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', 2);
            tc.verifyGreaterThan(g2.NumMembers, g1.NumMembers);
            % band members (2..end-1): peak wavenumber ratio between neighbours is 2^(-1/V)
            lq = linspace(0, max(lam), 20000)';
            for gv = {g1, g2}
                g = gv{1};  V = g.VoicesPerOctave;
                kp = zeros(1, g.NumMembers);
                H = graphfilters(g, 'Lambda', lq);
                for m = 1:g.NumMembers, [~, i] = max(H(:, m));  kp(m) = sqrt(lq(i)); end
                r = kp(2:end-2) ./ kp(3:end-1);
                tc.verifyEqual(r, 2^(-1/V) * ones(size(r)), 'RelTol', 2e-2, sprintf('V=%d', V));
            end
        end

        function edgeMembersCompleteTheFrame(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', 1);
            H = graphfilters(g, 'Lambda', [0; max(lam)]);
            tc.verifyEqual(H(1, 1), 1, 'AbsTol', 1e-12);                    % DC belongs to the low-pass
            tc.verifyEqual(H(2, end), 1, 'AbsTol', 1e-12);                  % lmax to the high-pass
            tc.verifyEqual(sum(H(1, 2:end).^2), 0, 'AbsTol', 1e-12);
        end

        function itIsParsevalOnASensorGraph(tc)
            % identity mass: summed over the members and the vertices, the coefficient
            % energy equals the field's energy, which is what the ingest's spatial bands
            % partition at the root of the sensor tree.
            arr = rheome.sensors.grid('Size', [5 6], 'Pitch', 1e-2);
            G = rheome.sensors.graph(arr, 'Faces', false);  M = rheome.sensors.modes(G);
            g = rheome.graphfilterbank(M.Lambda, 'Wavelet', 'logitersine', 'VoicesPerOctave', 2, 'Transform', M.T);
            rng(1);  x = randn(arr.nCh, 3);
            W = wt(g, x);
            tc.verifyEqual(sum(W.^2, 'all'), sum(x.^2, 'all'), 'RelTol', 1e-10);
            tc.verifyEqual(reshape(sum(W.^2, [1 3]), 1, []), sum(x.^2, 1), 'RelTol', 1e-10);   % per column too
        end

    end
end
