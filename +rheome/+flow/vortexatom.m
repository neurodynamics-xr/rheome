function out = vortexatom(name, varargin)
% FLOW.VORTEXATOM  A space-time vortex atom: spinning, non-stationary, inside one band.
%
%   out = rheome.flow.vortexatom('sub01')                       % an alpha vortex, 1.33 s long
%   out = rheome.flow.vortexatom(name, Band=[8 16], Chirality=-1)       % the mirror image
%   out = rheome.flow.vortexatom(name, Winding=-1, Phase=0)             % an antivortex; a source
%
% ⭐⭐ THE SPIN IS THE OSCILLATION, WHICH IS WHY THIS STAYS IN BAND. A carrier multiplied by an
% independent spin rate is TWO atoms at f0 +- f_spin, and those sidebands leave the octave as soon
% as the spin is quick: measured in-band energy falls to 2.5% at a 100 ms half turn and 0.1% at
% 50 ms (docs/2026-09-26-feature-table-design.md section 29). Here the pattern's rotation IS the
% band's oscillation -- one analytic atom, one peak -- so a fast spin is not a problem but the
% definition. ⭐ The two-atom construction fitted ~4 turns inside an alpha atom; this one fits 15,
% because the turn rate is fc rather than a separate dial bounded by the tile's width.
%
% THE CONSTRUCTION, a separable product of three factors:
%   z(v,t) = A(d_v) * exp(i*(m*theta_v + phase))  *  psi_fc(t - t0)
%   J(v,t) = Re(z)*e1(v) + Im(z)*e2(v)
% radial profile x angular harmonic x ANALYTIC temporal wavelet. The first two are
% rheome.flow.seedvortex; the third is one member of a timefilterbank, whose phase advance gives the
% spin and whose envelope gives the non-stationarity. Its bandwidth is the tile's by construction,
% so the atom cannot leave the band it was designed in.
%
% ⭐⭐ IT IS RANK TWO, AND NEVER MATERIALISED. Writing psi = a + i*b,
%       J(:,t) = a(t)*JA + b(t)*JB,      JB = the SAME field with Phase + pi/2
% so [3nV x nT] is never formed: the forward is (G*JA)*a + (G*JB)*b, two leadfield products for
% the whole record. ⭐ And the spatial quadrature partner is not a new object -- it is
% rheome.flow.seedvortex's existing Phase dial turned 90 degrees, the same dial that separates a source
% from a vortex. Use out.materialise() only for a figure.
%
% INPUTS
%   name    dataset, as rheome.load.bases
%   Vertex  seed (default seedvortex's, chosen against the gauge)   WavelengthMM (140)
%   Winding m, integer (default +1)                  Phase (pi/2 = vortex; 0 = source)
%   Band    [f_lo f_hi] Hz (default [8 16], alpha)   fc = the nearest bank member to sqrt(lo*hi)
%   Duration seconds (default 4)                     SampleRate Hz (default 600)
%   CentreSec  where the envelope peaks (default Duration/2)
%   Chirality  +1 (default) | -1   conjugates psi, reversing the sense of rotation
%   MomentNAm  TOTAL dipole moment in nAm; [] (default) leaves the peak-normalised field
%              ⚠ the total, NOT the per-vertex peak -- at 140 mm they differ by 652x
%   Hemi ("L")  Bases  Gauge  Check (true)
%
% OUTPUT (struct out)
%   .JA .JB  [3nV x 1]  the spatial quadrature pair, interleaved as rheome.flow.seedvortex
%   .a .b    [1 x nT]   real and imaginary parts of the temporal atom
%   .psi     [1 x nT]   the analytic temporal wavelet, envelope peak 1
%   .zs      [nV x 1]   the complex spatial factor        .fc  the member's centre frequency
%   .materialise  @() [3nV x nT], for figures only       .seed  rheome.flow.seedvortex's output
%   .check   peakHz, inbandFrac, supportSec, turnsInSupport, rotationHz, nStrong, charge,
%            momentPerUnitPeak, totalMomentNAm, peakPerVertexNAm
%
% ⭐ SCALE, AMPLITUDE AND RATE ARE INDEPENDENT. fc and the 1.26 s support are identical at 70, 140
%   and 267 mm, because the atom is a separable product -- changing the spatial factor cannot move
%   the temporal one. ⚠ But observability is NOT scale-independent: at a fixed 10 nAm TOTAL moment
%   the sensor peak falls 23.37 -> 17.73 -> 13.48 fT as the wavelength grows 70 -> 140 -> 267 mm,
%   so a coarse vortex is penalised twice, needing more total moment AND yielding less field per
%   unit of it. The |J| peak sits at lambda/2pi (11.0, 22.3, 42.5 mm measured).
%
% ⚠ peakHz is measured on the SOURCE field, per vertex. Whether the SENSORS see alpha is a
%   separate question and needs the leadfield -- rheome.flow.vortexatom does not answer it. (It does:
%   noiseless in-band 1.0000 at peak 11.25 Hz, and 0.6311 still peaking at 11.5 Hz with real
%   empty-room noise at SNR 3. See docs/2026-09-26-feature-table-design.md section 30.)
%
% ⚠⚠ THE RECORD MUST BE SEVERAL TIMES THE ATOM'S SUPPORT, and fc IS THE NEAREST MEMBER THAT
%   EXISTS, not the band's centre. Measured: asking for [4 8] on a 4 s record returns fc 6.727
%   rather than 5.657, because the bank has no 5.657 Hz member that short -- the grid truncates at
%   the low end -- and the support comes back 2.133 s, 53% of the record and truncated by the
%   periodic boundary. On 8 s the same request gives fc 5.657 and support 2.518 s, and the
%   constant-Q doubling against alpha's 1.262 s is then exact to 0.3% (1.995). `.check.truncated`
%   flags it and a warning fires. ⭐ Alpha's support is 1.26 s at every record length tried, so
%   4 s is ample for alpha and thin for anything below it.
%
% See also: rheome.flow.seedvortex, rheome.timefilterbank, rheome.selection.wavelettile, rheome.forward.leadfield,
%           rheome.spectral.aperiodic, rheome.operators.gauge
%
% Author: Diellor Basha, 2026

    p = inputParser;
    p.addParameter('Vertex', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('WavelengthMM', 140, @isscalar);
    p.addParameter('Winding', 1, @(x) isscalar(x) && x == fix(x));
    p.addParameter('Phase', pi/2, @isscalar);
    p.addParameter('Band', [8 16], @(x) numel(x)==2 && x(2)>x(1));
    p.addParameter('Duration', 4, @isscalar);
    p.addParameter('SampleRate', 600, @isscalar);
    p.addParameter('CentreSec', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('Chirality', 1, @(x) isscalar(x) && abs(x)==1);
    p.addParameter('MomentNAm', [], @(x) isempty(x) || (isscalar(x) && x > 0));
    p.addParameter('Hemi', "L");
    p.addParameter('Bases', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Gauge', [], @(x) isempty(x) || isstruct(x));
    p.addParameter('Check', true, @islogical);
    p.parse(varargin{:});
    o = p.Results;

    fs = o.SampleRate;
    nT = round(o.Duration * fs);
    if isempty(o.CentreSec), o.CentreSec = o.Duration/2; end

    %% 1. the spatial factors -- one seedvortex call, lifted twice
    % ⚠ the two lifts MUST share a normaliser. |1i*z| = |z| pointwise so the max is the same, but
    %   computing it once makes that explicit rather than incidental.
    sv = rheome.flow.seedvortex(name, Vertex=o.Vertex, WavelengthMM=o.WavelengthMM, ...
                         Winding=o.Winding, Phase=o.Phase, Hemi=o.Hemi, ...
                         Bases=o.Bases, Gauge=o.Gauge, Check=o.Check);
    z = sv.z;
    g = o.Gauge;
    if isempty(g)
        B = o.Bases;  if isempty(B), B = rheome.load.bases(name); end
        H = B.(char(o.Hemi));
        g = rheome.operators.gauge(H.S.Vertices, double(H.S.Faces), Method="diffusion");
    end
    % ⚠⚠ THE NORMALISER IS PEAK-BASED, SO IT IS NOT AN AMPLITUDE. max|J| = 1 means the TOTAL
    %   dipole moment depends on the scale: measured 147.9 / 652.6 / 2735.2 at 70 / 140 / 267 mm,
    %   roughly four-fold per octave because it follows the area. Comparing two scales at equal
    %   PEAK compares sources 18x apart in total moment. MomentNAm removes the trap by fixing the
    %   total instead, which is the quantity a physiological figure refers to.
    nrm = max(vecnorm(real(z).*g.e1 + imag(z).*g.e2, 2, 2));
    JA = reshape((  real(z).*g.e1 + imag(z).*g.e2 )', [], 1) / nrm;   % lift(z)
    JB = reshape(( -imag(z).*g.e1 + real(z).*g.e2 )', [], 1) / nrm;   % lift(1i*z) = Phase + pi/2

    momentPerUnitPeak = sum(vecnorm(reshape(JA,3,[])', 2, 2));
    if ~isempty(o.MomentNAm)
        amp = o.MomentNAm*1e-9 / momentPerUnitPeak;      % scale BOTH, or the circle turns elliptic
        JA = JA*amp;  JB = JB*amp;
    end

    %% 2. the temporal factor -- one analytic member of the constant-Q bank
    tfb = rheome.timefilterbank(nT, 'SamplingFrequency', fs);
    fcAll = centerFrequencies(tfb);
    f0 = sqrt(o.Band(1)*o.Band(2));
    [~, mIdx] = min(abs(fcAll(:) - f0));
    fc = fcAll(mIdx);
    [Hm, ~] = freqz(tfb);
    % ⚠ one-sided only: leaving the negative half at zero is what makes psi ANALYTIC, and an
    %   analytic atom is what carries a SIGNED rotation sense. Mirroring it would give a real
    %   wavelet whose two chiralities are indistinguishable.
    W = zeros(1, nT);
    nB = size(Hm, 2);
    W(1:nB) = Hm(mIdx, :);
    psi = ifft(W);
    psi = circshift(psi, round(o.CentreSec*fs));
    psi = psi / max(abs(psi));
    if o.Chirality < 0, psi = conj(psi); end
    a = real(psi);  b = imag(psi);

    out = struct('JA', JA, 'JB', JB, 'a', a, 'b', b, 'psi', psi, 'zs', z, ...
                 'fc', fc, 'band', o.Band, 'member', mIdx, 'nT', nT, 'fs', fs, ...
                 'winding', o.Winding, 'phase', o.Phase, 'chirality', o.Chirality, ...
                 'vertex', sv.vertex, 'seed', sv, 'check', struct());
    out.materialise = @() JA*a + JB*b;

    %% 3. verify, on the source field, that it is what it claims
    env = abs(psi);
    e2c = env.^2;  [~, ord] = sort(e2c, 'descend');
    keep = false(1, nT);  keep(ord(1:find(cumsum(e2c(ord)) >= 0.999*sum(e2c), 1))) = true;
    out.check.supportSec     = sum(keep)/fs;
    out.check.truncated      = out.check.supportSec > 0.5*o.Duration;
    if out.check.truncated
        warning('flow:vortexatom:truncated', ...
            ['The atom''s support is %.3f s in a %.3f s record (%.0f%%), so the periodic ' ...
             'boundary truncates it. Lengthen Duration to several times the support.'], ...
            out.check.supportSec, o.Duration, 100*out.check.supportSec/o.Duration);
    end
    out.check.turnsInSupport = out.check.supportSec * fc / abs(o.Winding);
    % the spin: unwrapped phase advance of the temporal atom over its support
    ph = unwrap(angle(psi(keep)));
    out.check.rotationHz = o.Chirality * abs(median(diff(ph))) * fs / (2*pi) / abs(o.Winding);
    % peak frequency and in-band share, on the strongest vertices' own time series
    amp = abs(z);  strong = amp > 0.5*max(amp);
    out.check.nStrong = sum(strong);
    Xs = real(z(strong) * psi);                       % [nS x nT], the e1 component
    P  = abs(fft(Xs, [], 2)).^2;
    ff = (0:nT-1)*fs/nT;  k = 1:floor(nT/2);
    [~, ip] = max(P(:, k), [], 2);
    out.check.peakHz = median(ff(k(ip)));
    inb = ff(k) >= o.Band(1) & ff(k) <= o.Band(2);
    out.check.inbandFrac = median(sum(P(:, k(inb)), 2) ./ sum(P(:, k), 2));
    out.check.momentPerUnitPeak = momentPerUnitPeak;
    out.check.totalMomentNAm    = sum(vecnorm(reshape(JA,3,[])', 2, 2))*1e9;
    out.check.peakPerVertexNAm  = max(vecnorm(reshape(JA,3,[])', 2, 2))*1e9;
    out.check.charge = NaN;
    if isfield(sv.check, 'charge'), out.check.charge = sv.check.charge; end
end

% Author: Diellor Basha, 2026
