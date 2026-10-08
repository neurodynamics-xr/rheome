classdef tGfbAtoms < matlab.unittest.TestCase

    methods (Test)

        function atomsMatchFiltersLocalize(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            b  = struct('Phi',fx.Phi, 'Lambda',lam, 'Mass',fx.Mass, 'nV',fx.nV);
            f  = rheome.filters.frame('mexhat', 6, lam, 'warn', false);
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',6, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            seeds = [1; 42; 300];
            ref = rheome.filters.localize(b, seeds, rheome.filters.frame_gains(f, lam));   % [nV x P x M]
            A   = graphfilters(g, 'Type','vertex', 'Vertex', seeds);
            tc.verifyEqual(A, ref, 'AbsTol', 1e-10);
        end

        function singleSeedKeepsItsDimension(tc)
            % rheome.filters.localize collapses [nV x 1 x M] to [nV x M]; we deliberately do not,
            % so the shape does not depend on how many seeds were asked for.
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',6, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            A = graphfilters(g, 'Type','vertex', 'Vertex', 17);
            tc.verifySize(A, [fx.nV, 1, g.NumMembers]);
        end

        function atomIsLocalisedAtItsSeed(tc)
            % An atom is a wavelet: on a sphere the kernel is radially symmetric about
            % the seed, so its peak belongs there. Member 2 is the FINEST wavelet
            % (member 1 is the low-pass scaling function; t ascends after it).
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',6, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            seed = 100;
            A = graphfilters(g, 'Type','vertex', 'Vertex', seed);
            [~, ipk] = max(abs(A(:,1,2)));
            d = norm(fx.V(ipk,:) - fx.V(seed,:));
            tc.verifyLessThan(d, 0.02);                % within 20 mm on a 100 mm sphere
        end

        function vertexTypeRequiresAVertex(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            tc.verifyError(@() graphfilters(g, 'Type','vertex'), 'graphfilterbank:vertex');
        end

        function vertexTypeRequiresATransform(tc)
            g = rheome.graphfilterbank(100);
            tc.verifyError(@() graphfilters(g, 'Type','vertex', 'Vertex', 1), ...
                'graphfilterbank:noTransform');
        end

        function impulseMatchesTheVertexModeOfGraphfilters(tc)
            % impulse is the SAME operation graphfilters('Type','vertex') performs, named
            % as jointfilterbank names it. If these ever disagree there are two
            % implementations of one idea, which is the thing this method exists to avoid.
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',6, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            v = 137;
            A = graphfilters(g, 'Type','vertex', 'Vertex', v);
            for m = 1:g.NumMembers
                tc.verifyEqual(impulse(g, m, v), A(:, 1, m), 'AbsTol', 1e-12);
            end
        end

        function impulseReturnsOneColumnPerVertex(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Wavelet','mexhat', 'NumFilters',6, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            tc.verifySize(impulse(g, 2, 137), [fx.nV, 1]);
        end

        function impulseRefusesAVertexOutsideTheGraph(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            tc.verifyError(@() impulse(g, 1, fx.nV + 1), 'graphfilterbank:vertex');
            tc.verifyError(@() impulse(g, 1, 2.5),       'graphfilterbank:vertex');
        end

        function impulseRefusesAMemberOutsideTheBank(tc)
            fx = gfbFixture();  lam = fx.Lambda;
            g  = rheome.graphfilterbank(lam, 'NumFilters',6, ...
                                 'Transform', rheome.graphtransform.eigen(fx.Phi, fx.Mass, lam));
            tc.verifyError(@() impulse(g, g.NumMembers + 1, 1), 'graphfilterbank:member');
        end

        function impulseRequiresATransform(tc)
            % A design-stage bank has no way to reach the vertex domain at all.
            g = rheome.graphfilterbank(100);
            tc.verifyError(@() impulse(g, 1, 1), 'graphfilterbank:noTransform');
        end

        function theChebyshevRouteProducesAnAtomAtAll(tc)
            % ⚠ THE CHEBYSHEV TRANSFORM HAS IDENTITY forward/inverse. Its filtering happens
            % through .filter, by polynomial recursion in L. Anything that reaches for
            % forward/inverse gets the delta back unchanged and multiplies it by a gain
            % vector, which is silently WRONG rather than an error.
            [g, ~, v] = tGfbAtoms.senbank(60);
            u = impulse(g, 4, v);
            tc.verifyGreaterThan(sum(u ~= 0), 1, 'a Chebyshev atom must not be all zeros');
        end

        function theTwoRoutesAgreeOnTheSameAtom(tc)
            % The claim the whole two-route design rests on: same operator, same member,
            % same atom, whichever way you reach the vertex domain.
            [gc, ge, v] = tGfbAtoms.senbank(60);
            uc = impulse(gc, 4, v);  ue = impulse(ge, 4, v);
            tc.verifyEqual(norm(uc - ue)/norm(ue), 0, 'AbsTol', 0.02);
        end

        function aPolynomialFilterTouchesExactlyItsOwnDegreeInHops(tc)
            % g(L) of degree K is a polynomial in L, and L only connects neighbours, so the
            % atom is EXACTLY zero beyond K hops. That is compact support, and it is the one
            % real difference between the two routes: the exact eigenbasis has none.
            [~, ge, v, G] = tGfbAtoms.senbank(60);
            hop = tGfbAtoms.hops(G, v);
            for K = [3 6]
                T  = rheome.graphtransform.chebyshev(G.L, max(G.Lambda_), 'Order', K);
                gK = rheome.graphfilterbank(G.Lambda_, 'Wavelet','itersine', 'NumFilters',6, ...
                                     'Transform', T);
                u  = impulse(gK, 4, v);
                tc.verifyEqual(max(hop(u ~= 0)), K, ...
                    sprintf('order %d must reach exactly %d hops', K, K));
                tc.verifyTrue(all(u(hop > K) == 0), 'beyond K hops must be exactly zero');
            end
            % the exact route reaches the whole graph, with no zero anywhere
            ue = impulse(ge, 4, v);
            tc.verifyEqual(sum(ue == 0), 0);
        end

    end

    methods (Static)
        function [gc, ge, v, G] = senbank(order)
            % A small array whose FULL spectrum is available, so lambda_max really does
            % bound the operator -- the condition rheome.graphtransform.chebyshev insists on.
            persistent P
            if isempty(P)
                w = warning('off', 'sensors:calibrate:narrowFitWindow');
                cl = onCleanup(@() warning(w));
                arr = rheome.sensors.grid('Size', [12 12], 'Pitch', 10e-3);
                Gs  = rheome.sensors.graph(arr);
                Ms  = rheome.sensors.modes(Gs);
                Gs.Lambda_ = Ms.Lambda;  Gs.Phi_ = Ms.Phi;
                c = mean(Gs.Vertices, 1);
                [~, iv] = min(vecnorm(Gs.Vertices - c, 2, 2));
                P = struct('G', Gs, 'v', iv);
            end
            G = P.G;  v = P.v;
            gc = rheome.graphfilterbank(G.Lambda_, 'Wavelet','itersine', 'NumFilters',6, ...
                     'Transform', rheome.graphtransform.chebyshev(G.L, max(G.Lambda_), 'Order', order));
            ge = rheome.graphfilterbank(G.Lambda_, 'Wavelet','itersine', 'NumFilters',6, ...
                     'Transform', rheome.graphtransform.eigen(G.Phi_, speye(G.nV), G.Lambda_));
        end

        function hop = hops(G, v)
            A = double(G.W ~= 0);  A(1:G.nV+1:end) = 0;
            hop = inf(G.nV, 1);  hop(v) = 0;  front = v;  h = 0;
            while ~isempty(front) && h < G.nV
                h = h + 1;
                nxt = find(any(A(:, front) ~= 0, 2) & ~isfinite(hop));
                hop(nxt) = h;  front = nxt;
            end
        end
    end
end

% Author: Diellor Basha, 2026
