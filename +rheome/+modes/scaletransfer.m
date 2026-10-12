function Tr = scaletransfer(S, win, opts)
% MODES.SCALETRANSFER  Per participant: how much of the energy planted at graph scale s the instrument returns at s.
%
%   Tr = rheome.modes.scaletransfer(S, win)       % win: table depth, lo_mm, hi_mm, scale_mm (rheome.modes.jointfields)
%   Tr = rheome.modes.scaletransfer(S, win, Edge=0.5)
%
% Protocol section 14, J4 (resolution per measurement, never a cut), extending the modulation transfer of
% rheome.inverse.resolution (.mtf: the diagonal of the resolution matrix R in the normally oriented mode family,
% B_kk = <n phi_k, R n phi_k>, the amplitude kernel's Wiener gain per mode) from modes to the graph scales of
% the joint field. A random normal current is planted in one hemisphere with its energy in scale window s only
% (Laplace-Beltrami coefficients white, weighted by g_s, the window's members of the tight logitersine bank that
% rheome.modes.operator uses), seen through the participant's gain and unconstrained minimum-norm kernel, and
% read with the same bank. Exact in expectation, no draws (B = Phi' M Kq Gq Phi):
%
%   T_abs(s)  E<P_s q, P_s q_hat> / E||P_s q||^2 = sum_k g_s^4 B_kk / sum_k g_s^4: the share of the planted
%             scale-s pattern returned, coherently, at s -- the MTF averaged over the scale with the bank's
%             weights. Absolute (1 = perfect), small for normally oriented patterns (a smooth normal pattern
%             flips across every gyrus and largely cancels at the sensors: rheome.inverse.resolution)
%   T(s)      T_abs / its peak over every member of the bank (the finer members included): the curve the
%             scale edge is read on, as .mtfHalf reads half of the MTF's peak. Each atom's T-bar uses T
%   edge_mm   where T falls through Edge (0.5) going finer from its peak over the depths, interpolated in log
%             scale; NaN with edge_flag "beyond_finest" (T >= Edge from the peak to the finest depth) or
%             "below_edge" (no depth reaches Edge). The coarsest depth sits below its neighbour (one OMEGA head:
%             T 0.14 at 400 mm, 1.00 at 270 mm; the near-constant mode cancels at the sensors) and is reported,
%             not used for the edge. Same head: edge 114 mm (T 0.55 at 130 mm, 0.42 at 92 mm). ⚠ Like .mtfHalf this
%             crossing can jump on a bumpy curve: centroid_mm (the T-weighted centroid in log scale over the
%             depths, as .mtfCentroid) is the stable summary
%   share(s)  purity: of the estimate's whole energy (both hemispheres, every member), the share that falls
%             in window s of the planted hemisphere, over the same share for a perfect instrument (B = I).
%             ⚠ Measured on one OMEGA head: 0.03-0.14 at depths 0-5 -- whatever the planted scale, about
%             70 % of the minimum-norm normal current's energy lies in the members finer than 56 mm
%
% ⚠ Noiseless: the resolution part of the transfer; SNR is reported separately. ⚠ Energy outside the K
% eigenmodes of the hemisphere's basis is not seen by the bank and is left out.
%
% Measured on that head: T_abs peaks at 0.019 (270 mm).
%
% Tr.table (depth, scale_mm, T, T_abs, share, and per hemisphere T_abs_<h>, share_<h>)  Tr.member (wavelength_mm,
% T_abs per bank member, mean of hemispheres)  Tr.peak  Tr.edge_mm  Tr.edge_flag  Tr.centroid_mm
%
% See also: rheome.modes.jointfields, rheome.modes.operator, rheome.inverse.resolution
%
% Author: Diellor Basha, 2026

    arguments
        S (1,1) struct
        win table
        opts.Voices (1,1) double = 2
        opts.Hemis string = ["L" "R"]
        opts.Edge (1,1) double = 0.5
    end
    K = S.Res.ImagingKernel;  nr = double(S.Sf.VertNormals);  nr = nr ./ vecnorm(nr, 2, 2);
    nh = numel(opts.Hemis);  H = struct([]);
    for i = 1:nh
        Hh = S.B.(opts.Hemis(i));  gv = double(Hh.gv(:));  lbo = Hh.lbo;  n = nr(gv, :);
        Kq = n(:,1) .* K((gv-1)*3 + 1, :) + n(:,2) .* K((gv-1)*3 + 2, :) + n(:,3) .* K((gv-1)*3 + 3, :);
        Gq = n(:,1)' .* S.G(:, (gv-1)*3 + 1) + n(:,2)' .* S.G(:, (gv-1)*3 + 2) + n(:,3)' .* S.G(:, (gv-1)*3 + 3);
        lam = lbo.Lambda(:);
        gfb = rheome.graphfilterbank(lam, 'Wavelet', 'logitersine', 'VoicesPerOctave', opts.Voices);
        G2 = cell2mat(arrayfun(@(m) feval(gain(gfb, m), lam), 1:gfb.NumMembers, 'uni', 0)).^2;
        wl = 1e3 * wavelengths(gfb);
        H(i).B = lbo.Phi' * (lbo.Mass * Kq);  H(i).C = Gq * lbo.Phi;  H(i).G2 = G2;  H(i).wl = wl;
        H(i).in = arrayfun(@(s) wl >= win.lo_mm(s) & wl < win.hi_mm(s), 1:height(win), 'uni', 0);
    end
    ns = height(win);  Ta = nan(ns, nh);  sh = nan(ns, nh);  nm = size(H(1).G2, 2);  Tm = nan(nm, nh);
    for i = 1:nh
        BX = arrayfun(@(j) (H(j).B * H(i).C).^2, 1:nh, 'uni', 0);       % [K_j x K_i] squared transfer, i -> j
        Bd = sum(H(i).B .* H(i).C', 2);                                   % B_kk: the MTF of rheome.inverse.resolution
        if size(H(i).G2, 2) == nm, Tm(:, i) = (H(i).G2.^2)' * Bd ./ sum(H(i).G2.^2, 1)'; end
        for s = 1:ns
            gs2 = sum(H(i).G2(:, H(i).in{s}), 2);
            if ~any(gs2), continue; end
            Ta(s, i) = sum(gs2.^2 .* Bd) / sum(gs2.^2);
            tot = 0;
            for j = 1:nh
                Em = H(j).G2' * (BX{j} * gs2);  tot = tot + sum(Em);
                if j == i, ret = sum(Em(H(i).in{s})); end
            end
            Ei = H(i).G2' * gs2;  sh(s, i) = (ret / tot) / (sum(Ei(H(i).in{s})) / sum(Ei));
        end
    end
    pk = max(mean(Tm, 2, 'omitnan'));  Tabs = mean(Ta, 2, 'omitnan');  T = Tabs / pk;
    tab = table(win.depth, win.scale_mm, T, Tabs, mean(sh, 2, 'omitnan'), 'VariableNames', {'depth','scale_mm','T','T_abs','share'});
    for i = 1:nh, tab.("T_abs_" + opts.Hemis(i)) = Ta(:, i);  tab.("share_" + opts.Hemis(i)) = sh(:, i); end
    [~, o] = sort(tab.scale_mm, 'descend');  t = T(o);  sc = tab.scale_mm(o);
    [tp, ip] = max(t);  j = find(t(ip:end-1) >= opts.Edge & t(ip+1:end) < opts.Edge, 1) + ip - 1;
    if tp < opts.Edge, e = NaN;  flag = "below_edge";
    elseif isempty(j), e = NaN;  flag = "beyond_finest";
    else, e = exp(interp1([t(j) t(j+1)], log([sc(j) sc(j+1)]), opts.Edge));  flag = "interpolated";
    end
    ok = isfinite(t);  cen = exp(sum(t(ok) .* log(sc(ok))) / sum(t(ok)));
    Tr = struct('table', tab, 'member', table(H(1).wl(:), mean(Tm, 2, 'omitnan'), 'VariableNames', {'wavelength_mm','T_abs'}), ...
                'peak', pk, 'edge_mm', e, 'edge_flag', flag, 'centroid_mm', cen, 'edge', opts.Edge);
end

% Author: Diellor Basha, 2026
