function out = vortextrack(bank, D, varargin)
% FLOW.VORTEXTRACK  A time-resolved PARAMETRIC inverse: vortex parameters as functions of time.
%
%   out = rheome.flow.vortextrack(bank, D)
%   out = rheome.flow.vortextrack(bank, D, Cycles=[2 5 15], Stride=0.1)
%
% ⭐⭐ WHAT MAKES THIS AN INVERSE RATHER THAN REPEATED MATCHING. A distributed inverse solves for
% 3nV unknowns per sample and regularises the result into existence. This solves for SIX numbers per
% time point -- place, scale, duration, amplitude, rotation and fit -- by choosing the atom that best
% explains the record there. The model is the parameterisation, so the answer is a trajectory of
% physical quantities rather than a map that still has to be interpreted.
%
% ⚠⚠ Q IS ONLY IDENTIFIABLE WITH A LENGTH PENALTY. Selecting on raw energy returns the LONGEST
% atom in the list regardless of the truth -- 0 of 6 correct for planted Q of 2 and of 5 at SNR 3.
% `Penalty` 0.5, dividing by sqrt(support), gets 17 of 18. See the code note.
%
% ⭐ Q IS A FITTED PARAMETER, NOT A SETTING. Section 34.1 measured that no single Q describes the
% alpha envelope: short atoms give a fast white envelope (tau 0.04 s, chi 0.13) and long ones a slow
% structured one (tau 0.27, chi 1.33), while the recording needs 0.19 and 0.78 at once. Fitting Q per
% time point is the parametric form of the mixture that resolves it, and `out.cycles` is then a
% measurement -- how many cycles the rhythm sustained at that moment -- not a choice.
%
% ⭐ IT COSTS TWO CONVOLUTIONS AND TWO PRODUCTS PER Q. The spatial factors do not depend on the
% temporal one, so one bank serves every Q. Sliding the atom is a convolution, so
%     U_q = D * conj(a_q),  V_q = D * conj(b_q)         [nCh x nT], by FFT
%     c0  = PA'*U_q + PB'*V_q,  c1 = PB'*U_q - PA'*V_q  [nAtoms x nT]
% gives the phase-invariant energy of every atom at every instant. Nothing of size [nCh x nT] is
% formed per atom.
%
% ⚠ THE TRACK IS NOT A DETECTION. `energyFrac` says how much of the record the winning atom
%   explains, not whether a vortex is there: section 32 measured a standing field scoring 0.697
%   against a rotating one's 0.900. Read `icoh` alongside it, which separated the two at every
%   location tested, and treat a high energy with a near-zero icoh as a focal source that is not
%   turning.
%
% INPUTS
%   bank    rheome.flow.vortexbank output          D  [nCh x nT] sensor record, bank.chan order
%   Cycles  Q values to fit ([2 5 15])      Stride  seconds between estimates (0.05)
%   Fc      centre frequency to fit at, [] = the bank's own. The spatial dictionary is
%           band-independent, so one bank tracks every band.
%   Penalty selection divides energy by support^Penalty (0.5; see the note in the code)
%   Icoh    true (default) computes the rotation statistic at each estimate
%
% OUTPUT (struct out)
%   .t [1 x nE] seconds                    .table  one row per estimate
%   .vertex .scaleMM .cycles .energy .energyFrac .momentNAm .icoh   [nE x 1]
%   .energyAll [nAtoms x nCycles x nE] if KeepAll
%
% See also: rheome.flow.vortexbank, rheome.flow.vortexmatch, rheome.flow.vortexatom, rheome.flow.vortexspectrum
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Cycles', [2 5 15], @(x) isnumeric(x) && all(x > 0));
    p.addParameter('Stride', 0.05, @isscalar);
    p.addParameter('Icoh', true, @islogical);
    % ⭐⭐ THE LENGTH PENALTY IS NOT A TASTE PARAMETER, IT IS WHAT MAKES Q IDENTIFIABLE. Raw energy
    %   does not merely favour long atoms, it is PINNED to the longest one: over six draws at SNR 3
    %   it returned Q = 15 for planted Q of 2 and 5 every single time, 0/6 correct, because a longer
    %   atom overlaps more of the record and accumulates more energy whatever the truth. Dividing by
    %   sqrt(support) recovers 17 of 18, and dividing by support itself over-corrects to the
    %   shortest atom, 0/6 for planted 5 and 15. ⚠ Changing this means re-measuring that sweep.
    p.addParameter('Penalty', 0.5, @isscalar);
    % ⭐⭐ ONE BANK SERVES EVERY BAND. rheome.flow.seedvortex has no temporal factor, so PA/PB are purely
    %   spatial and do not depend on frequency at all. Only the kernel's centre moves, which makes a
    %   broadband parametric inverse cost two correlations per (band, Q) on a dictionary built once.
    p.addParameter('Fc', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('KeepAll', false, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    [nCh, nT] = size(D);
    assert(nCh == size(bank.PA,1), ...
        'rheome.flow.vortextrack: D has %d channels, the bank has %d.', nCh, size(bank.PA,1));
    fs = bank.fs;  nQ = numel(o.Cycles);
    fc = o.Fc;  if isempty(fc), fc = bank.fc; end
    assert(fc > 0 && fc < fs/2, 'rheome.flow.vortextrack: Fc %g is not inside (0, %g).', fc, fs/2);
    N  = size(bank.PA,2);

    tt = (-fix(nT/2):fix((nT-1)/2))/fs;
    step = max(1, round(o.Stride*fs));
    idx  = 1:step:nT;
    nE   = numel(idx);

    E = zeros(N, nQ, nE);
    for q = 1:nQ
        sd  = o.Cycles(q)/(2*pi*fc);
        psi = exp(2i*pi*fc*tt) .* exp(-tt.^2/(2*sd^2));
        psi = ifftshift(psi);  psi = psi/norm(psi);
        a = real(psi);  b = imag(psi);
        % ⚠ correlation, not convolution: conj and reverse, or every estimate lands mirrored in time
        Ka = conj(fft(a));  Kb = conj(fft(b));
        Fd = fft(D, [], 2);
        U  = real(ifft(Fd .* Ka, [], 2));
        V  = real(ifft(Fd .* Kb, [], 2));
        aa = a*a';  bb = b*b';  ab = a*b';
        n0 = aa*sum(bank.PA.^2,1)' + bb*sum(bank.PB.^2,1)' + 2*ab*sum(bank.PA.*bank.PB,1)';
        n1 = aa*sum(bank.PB.^2,1)' + bb*sum(bank.PA.^2,1)' - 2*ab*sum(bank.PA.*bank.PB,1)';
        c0 = bank.PA'*U(:,idx) + bank.PB'*V(:,idx);
        c1 = bank.PB'*U(:,idx) - bank.PA'*V(:,idx);
        E(:,q,:) = reshape(c0.^2./n0 + c1.^2./n1, N, 1, nE);
    end

    wq = reshape(o.Cycles(:).^(-o.Penalty), 1, nQ, 1);
    [~, lin] = max(reshape(E.*wq, N*nQ, nE), [], 1);
    [ia, iq] = ind2sub([N nQ], lin);
    emax = zeros(1, nE);
    for e = 1:nE, emax(e) = E(ia(e), iq(e), e); end     % report the RAW energy of the winner

    % local data energy in a window of the winning atom's own length, for a comparable fraction
    frac = zeros(1,nE);  mom = zeros(1,nE);
    for e = 1:nE
        half = round(o.Cycles(iq(e))/fc*fs);
        lo = max(1, idx(e)-half);  hi = min(nT, idx(e)+half);
        den = sum(sum(D(:,lo:hi).^2));
        frac(e) = emax(e)/max(den, realmin);
        % ⭐ the moment follows from the fit, since the forward is linear in it
        k = ia(e);
        gain = sqrt(sum(bank.PA(:,k).^2) + sum(bank.PB(:,k).^2));
        mom(e) = sqrt(max(emax(e),0))/max(gain,realmin)*1e9;
    end

    ic = nan(1,nE);
    if o.Icoh
        for e = 1:nE
            k = ia(e);
            half = max(4, round(1.5*o.Cycles(iq(e))/fc*fs));
            lo = max(1, idx(e)-half);  hi = min(nT, idx(e)+half);
            if hi-lo < 8, continue; end
            pA = bank.PA(:,k);  pB = bank.PB(:,k);
            xA = (pA'*D(:,lo:hi))/max(pA'*pA, realmin);
            xB = (pB'*D(:,lo:hi))/max(pB'*pB, realmin);
            zA = hilbert(xA(:));  zB = hilbert(xB(:));
            ic(e) = imag(mean(zA.*conj(zB))) / ...
                    max(sqrt(mean(abs(zA).^2)*mean(abs(zB).^2)), realmin);
        end
    end

    out.fc         = fc;
    out.t          = (idx-1)/fs;
    out.vertex     = bank.vertex(ia);
    out.scaleMM    = bank.scaleMM(ia);
    out.cycles     = o.Cycles(iq)';
    out.energy     = emax(:);
    out.energyFrac = frac(:);
    out.momentNAm  = mom(:);
    out.icoh       = ic(:);
    out.table = table(out.t(:), out.vertex(:), out.scaleMM(:), out.cycles(:), ...
                      out.energyFrac(:), out.momentNAm(:), out.icoh(:), ...
        'VariableNames', {'t','vertex','scaleMM','cycles','energyFrac','momentNAm','icoh'});
    if o.KeepAll, out.energyAll = E; end
end

% Author: Diellor Basha, 2026
