function out = vortexmatch(bank, D, varargin)
% FLOW.VORTEXMATCH  Match a sensor record against a vortex dictionary, energy AND rotation.
%
%   out = rheome.flow.vortexmatch(bank, D)              % D [nCh x nT], channels as bank.chan
%   out = rheome.flow.vortexmatch(bank, D, Top=5)
%
% ⭐ THE MATCH IS PHASE-INVARIANT AND COSTS TWO PRODUCTS. An atom's family over temporal phase
% spans {M0, M1} with M0 = pA*a + pB*b and M1 = pB*a - pA*b, so the best energy over all phases is
% <D,M0>^2/|M0|^2 + <D,M1>^2/|M1|^2, and both inner products come from u = D*a' and v = D*b'. No
% atom is ever synthesised as [nCh x nT].
%
% ⚠⚠ ENERGY LOCATES; IT DOES NOT ESTABLISH ROTATION. At SNR 3 over 40 locations, best match energy
% as a fraction of data energy is 0.900 for a rotating vortex (a very tight 0.899-0.901) but 0.697
% for a STANDING one, ranging 0.232 to 0.885. So it separates structure from empty-room noise by a
% factor near 9000, and a rotating field from a standing one by a per-location ratio of only 1.29
% in the median, under 1.5 at 28 of the 40 locations. ⚠ And the spread is wide enough that one
% seed settles nothing: the same comparison gives 1.02 at one location and 3.87 at another.
% ⭐ .icoh is the discriminator that does reach rotation: the imaginary coherency between the
% record's projections onto pA and pB, blind to the instantaneous real-valued mixing that leakage
% introduces. It separated the two at every location tested -- median 0.801 rotating against 0.002
% standing -- though the worst case is thinner than that headline, 0.080 against 0.023.
% Read .energy to say WHERE, and .icoh to say WHETHER IT TURNS, and never .energy as both.
%
% ⚠⚠ .icoh IS SIGNED AND THE SIGN IS ONLY THE HANDEDNESS WHEN THE LOCATION IS ALREADY PINNED.
%   Measured at SNR 3, 4 sites x 2 chiralities x 8 draws: the sign is right 64/64 when the true
%   vertex is ON the atom grid and 40/64 when it is not. ⚠ The off-grid failures are DETERMINISTIC
%   rather than noisy -- 0 of 8 at chirality -1 at three of the four sites, 8 of 8 at +1 at all
%   four -- so a coarse bank does not degrade toward chance, it reports one direction regardless,
%   which at the population level manufactures a preferred sense of rotation out of nothing.
%   The reason: .icoh is the imaginary coherency of pA'*D against pB'*D, and for D = u*a +/- w*b
%   the cross terms pA'*w and pB'*u are small only when the atom sits at the true vertex. Off-grid
%   its quadrature frame is rotated and the frame offset sets the sign instead. ⭐ So read the sign
%   after localising -- from a bank refined around the match -- never from the coarse top row.
%   Reporting |icoh| is the safe default and throws away the one quantity a phase-invariant energy
%   cannot recover. See docs/2026-09-26-feature-table-design.md section 54.
%
% INPUTS   bank  rheome.flow.vortexbank output    D  [nCh x nT]    Top  how many atoms to return (5)
%
% OUTPUT (struct out)
%   .table   Top rows: rank, vertex, scaleMM, energy, energyFrac, icoh
%   .best    the top row as a struct        .energy [N x 1] every atom's energy
%   .icoh    [N x 1] for the returned rows only (NaN elsewhere -- the Hilbert costs more)
%
% See also: rheome.flow.vortexbank, rheome.flow.vortexatom, rheome.spectral.carrier
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Top', 5, @(x) isscalar(x) && x >= 1);
    p.parse(varargin{:});
    o = p.Results;

    assert(size(D,1) == size(bank.PA,1), ...
        'rheome.flow.vortexmatch: D has %d channels, the bank has %d.', size(D,1), size(bank.PA,1));
    assert(size(D,2) == bank.nT, ...
        'rheome.flow.vortexmatch: D has %d samples, the bank was built for %d.', size(D,2), bank.nT);

    u = D*bank.a';  v = D*bank.b';
    c0 = bank.PA'*u + bank.PB'*v;
    c1 = bank.PB'*u - bank.PA'*v;
    E  = c0.^2./bank.n0 + c1.^2./bank.n1;

    [~, ord] = sort(E, 'descend');
    top = ord(1:min(round(o.Top), numel(ord)));
    dE  = norm(D, 'fro')^2;

    ic = nan(numel(E), 1);
    for k = top(:)'
        pA = bank.PA(:,k);  pB = bank.PB(:,k);
        xA = (pA'*D) / max(pA'*pA, realmin);        % the record's projection on each map
        xB = (pB'*D) / max(pB'*pB, realmin);
        zA = hilbert(xA(:));  zB = hilbert(xB(:));
        % ⭐ imaginary coherency: blind to instantaneous mixing, hence to leakage
        ic(k) = imag(mean(zA.*conj(zB))) / max(sqrt(mean(abs(zA).^2)*mean(abs(zB).^2)), realmin);
    end

    out.energy = E;
    out.icoh   = ic;
    out.table  = table((1:numel(top))', bank.vertex(top), bank.scaleMM(top), E(top), ...
                       E(top)/max(dE, realmin), ic(top), ...
        'VariableNames', {'rank','vertex','scaleMM','energy','energyFrac','icoh'});
    out.best   = table2struct(out.table(1,:));
end

% Author: Diellor Basha, 2026
