function out = vortexpursuit(bank, D, varargin)
% FLOW.VORTEXPURSUIT  Matching pursuit over the vortex dictionary: the parametric inverse.
%
%   out = rheome.flow.vortexpursuit(bank, D, NAtoms=40)
%   out = rheome.flow.vortexpursuit(bank, D, NAtoms=40, Source=true, Name="sub01")
%
% Greedy: fit the best atom over all (place, scale, duration, time), subtract its least-squares
% projection, repeat. The result is a list of atoms with parameters attached and a reconstruction
% built from them, which is the parametric counterpart of a distributed inverse's current map.
%
% ⭐ THE COMPARISON WITH A DISTRIBUTED INVERSE IS ABOUT PARAMETER COUNT. rheome.inverse.mne solves 3nV
% unknowns per sample under a norm penalty; this solves 6 per atom for a handful of atoms. A
% sensor-space goodness-of-fit therefore CANNOT compare them -- MNE fits the sensors almost exactly
% by construction and would win every time without that meaning anything. ⭐⭐ The honest test is
% prediction of data the fit never saw: hold out channels, fit both on the rest, and compare the
% predicted field at the held-out sensors. `Channels` does exactly that.
%
% ⚠ IT IS A BURST MODEL. Measured on 30 s of real alpha, 40 atoms explain 24.6% of band variance
% but the estimated-to-observed GFP ratio is 0.43 in the loudest quartile and 0.04 in the quietest,
% and 88 of 270 channels come out WORSE than predicting zero, the worst at -133%. Matching pursuit
% minimises TOTAL residual energy, which does not imply per-channel improvement. See
% docs/2026-09-26-feature-table-design.md section 36.1.
%
% INPUTS
%   bank     rheome.flow.vortexbank output        D  [nCh x nT]
%   NAtoms   how many to fit (40)          Cycles  Q values ([2 5 15])   Fc  centre ([] = bank.fc)
%   Stride   seconds between candidate positions (0.5)
%   Channels rows of D to FIT on; the rest are left for prediction ([] = all)
%   Source   true also returns the [3nV x nT] source estimate (false; needs Name)
%   Name     dataset, required when Source is true
%   Bases, Gauge, Hemi   passed through when regenerating source fields
%
% OUTPUT (struct out)
%   .atoms   table: iter, t, vertex, scaleMM, cycles, c0, c1, residFrac
%   .Bhat    [nCh x nT] sensor reconstruction (all channels, including any held out)
%   .resid   [1 x NAtoms] residual energy fraction after each atom
%   .J       [3nV x nT] source estimate, when Source is true
%
% See also: rheome.flow.vortexbank, rheome.flow.vortextrack, rheome.flow.vortexmatch, rheome.inverse.mne
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('NAtoms', 40, @(x) isscalar(x) && x >= 1);
    p.addParameter('Cycles', [2 5 15], @isnumeric);
    p.addParameter('Fc', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Stride', 0.5, @isscalar);
    p.addParameter('Channels', [], @isnumeric);
    p.addParameter('Source', false, @islogical);
    p.addParameter('Name', "", @(x) isstring(x) || ischar(x));
    p.addParameter('Hemi', "L");
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Gauge', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Verbose', false, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    [nCh, nT] = size(D);
    assert(nCh == size(bank.PA,1), ...
        'rheome.flow.vortexpursuit: D has %d channels, the bank has %d.', nCh, size(bank.PA,1));
    ch = o.Channels;  if isempty(ch), ch = 1:nCh; end
    ch = ch(:)';
    fs = bank.fs;
    fc = o.Fc;  if isempty(fc), fc = bank.fc; end
    K  = round(o.NAtoms);

    % a bank restricted to the fitting channels, so held-out rows never enter the fit
    sub = bank;  sub.PA = bank.PA(ch,:);  sub.PB = bank.PB(ch,:);

    tt = (-fix(nT/2):fix((nT-1)/2))/fs;
    R  = D(ch,:);                                   % the residual, on the fitted channels only
    e0 = sum(R(:).^2);
    Bhat = zeros(nCh, nT);
    rec  = zeros(1, K);
    it_ = zeros(K,1); tA = it_; vA = it_; sA = it_; qA = it_; c0A = it_; c1A = it_;

    for k = 1:K
        tr = rheome.flow.vortextrack(sub, R, 'Cycles', o.Cycles, 'Stride', o.Stride, ...
                              'Fc', fc, 'Icoh', false);
        [~, be] = max(tr.energy);
        ai = find(bank.vertex == tr.vertex(be) & bank.scaleMM == tr.scaleMM(be), 1);
        ncy = tr.cycles(be);  ti = round(tr.t(be)*fs) + 1;

        psi = i_atom(tt, fc, ncy, ti);
        a = real(psi);  b = imag(psi);
        % least-squares coefficients on the FITTED channels, applied to all of them
        m0f = sub.PA(:,ai)*a + sub.PB(:,ai)*b;
        m1f = sub.PB(:,ai)*a - sub.PA(:,ai)*b;
        c0 = sum(sum(R.*m0f))/max(sum(m0f(:).^2), realmin);
        c1 = sum(sum(R.*m1f))/max(sum(m1f(:).^2), realmin);
        R  = R - c0*m0f - c1*m1f;
        Bhat = Bhat + c0*(bank.PA(:,ai)*a + bank.PB(:,ai)*b) ...
                    + c1*(bank.PB(:,ai)*a - bank.PA(:,ai)*b);

        rec(k) = sum(R(:).^2)/max(e0, realmin);
        it_(k)=k; tA(k)=tr.t(be); vA(k)=bank.vertex(ai); sA(k)=bank.scaleMM(ai);
        qA(k)=ncy; c0A(k)=c0; c1A(k)=c1;
        if o.Verbose && mod(k,10)==0
            fprintf('  pursuit %d/%d, residual %.3f\n', k, K, rec(k));
        end
    end

    out.atoms = table(it_, tA, vA, sA, qA, c0A, c1A, rec(:), ...
        'VariableNames', {'iter','t','vertex','scaleMM','cycles','c0','c1','residFrac'});
    out.Bhat  = Bhat;
    out.resid = rec;
    out.fc    = fc;

    if o.Source
        assert(strlength(string(o.Name)) > 0, ...
            'rheome.flow.vortexpursuit: Source=true needs Name, to regenerate the seed fields.');
        Bs = o.Bases;  if isempty(Bs), Bs = rheome.load.bases(char(o.Name)); end
        g  = o.Gauge;
        Hm = Bs.(char(o.Hemi));
        if isempty(g)
            g = rheome.operators.gauge(Hm.S.Vertices, double(Hm.S.Faces), Method="diffusion");
        end
        nV = size(Hm.S.Vertices,1);
        out.J = zeros(3*nV, nT);
        cache = containers.Map('KeyType','char','ValueType','any');
        for k = 1:K
            key = sprintf('%d_%g', vA(k), sA(k));
            if isKey(cache, key)
                JJ = cache(key);
            else
                sv = rheome.flow.seedvortex(char(o.Name), Vertex=vA(k), WavelengthMM=sA(k), ...
                         Hemi=o.Hemi, Bases=Bs, Gauge=g, Check=false);
                z  = sv.z;
                nr = max(vecnorm(real(z).*g.e1 + imag(z).*g.e2, 2, 2));
                JJ = {reshape((  real(z).*g.e1 + imag(z).*g.e2 )', [], 1)/nr, ...
                      reshape(( -imag(z).*g.e1 + real(z).*g.e2 )', [], 1)/nr};
                cache(key) = JJ;
            end
            psi = i_atom(tt, fc, qA(k), round(tA(k)*fs)+1);
            a = real(psi);  b = imag(psi);
            out.J = out.J + c0A(k)*(JJ{1}*a + JJ{2}*b) + c1A(k)*(JJ{2}*a - JJ{1}*b);
        end
    end
end

function psi = i_atom(tt, fc, ncyc, ti)
    sd  = ncyc/(2*pi*fc);
    psi = ifftshift(exp(2i*pi*fc*tt).*exp(-tt.^2/(2*sd^2)));
    psi = psi/norm(psi);
    psi = circshift(psi, ti-1);
end

% Author: Diellor Basha, 2026
