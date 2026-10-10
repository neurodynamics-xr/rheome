function cfg = config(opts)
% INGEST.CONFIG  Parameters of a constant-Q tile summary, in canonical order, with a hash.
%
%   cfg = rheome.ingest.config()
%   cfg = rheome.ingest.config(FrameFloor=0.25, Wavelet="morse", VoicesPerOctave=[], ...
%                       FrequencyLimits=[], Precision="single", MaxBytes=2.5e9, ...
%                       Paged="auto", PageOverhead=0.25, MinPageLength=16, MaxSupport=Inf)
%
% The struct is the whole configuration of a store: same recording + same cfg + same
% release give byte-identical values. Its fields are assigned in one fixed order so that
% the serialisation, and therefore .hash, does not depend on the order the options were
% given in. The hash is the first 8 hex digits of MD5 over jsonencode of the ordered
% fields; it goes into the store's meta and is what a rebuild is compared against.
%
% ⚠ FrequencyLimits = [] means THE BANK'S OWN RANGE (cwtfreqbounds for this record):
% both ends are where the wavelet fits, under Nyquist and inside the record. Setting
% limits truncates the spectrum by hand, which the design avoids (design §3.2).
%
% ⚠ FrameFloor*fs must be an integer; that is checked where fs is known (rheome.ingest.grid).
%
% INPUTS (name-value):
%   FrameFloor       level-0 frame length in seconds (0.25)
%   Wavelet          "morse" | "amor" | "bump", passed to cwtfilterbank ("morse")
%   VoicesPerOctave  scales per octave, and therefore scales per band ([] = 4 for the
%                    frame bank, 10 for Morse; see the note in the code)
%   FrequencyLimits  [] for the bank's own range, or [fLo fHi] Hz
%   Precision        "single" | "double" for the CWT buffer; sums are always double
%   MaxBytes         byte guard on the per-channel CWT WORKING SET, 5x the coefficient
%                    buffer as measured (2.5e9): whole-record at 600 Hz/600 s is 2.2e9, the
%                    2400 Hz native store 9.9e9 and therefore paged; its slow-band
%                    whole-record job (2.9e9) is split in two. Peak RSS is this plus
%                    MATLAB's ~2 GB baseline plus the channel chunk.
%   Paged            "auto" (page when the whole-record buffer exceeds MaxBytes) |
%                    "always" | "never"
%   PageOverhead     halo budget: a band's page core is at least 2*tSupport/PageOverhead (0.25)
%   MinPageLength    shortest page core in seconds (16); rounded up to a grid level
%   Bank             "frame" (timefilterbank: designed tight constant-Q frame, sub-band
%                    evaluation; the default) | "morse" (cwtfilterbank at the full rate;
%                    the reference, 7x slower)
%   Anchor           frame bank only: absolute grid anchor in Hz (1); member k at Anchor*2^(k/V),
%                    bands are the integer octaves of the anchor
%   Oversample       frame bank only: sub-band samples per bin (2)
%   Space            also store the wavelength axis: a tight log-itersine frame on the
%                    sensor graph's spectrum applied to the field at every sample (true;
%                    needs sensor positions, otherwise skipped)
%   SpaceVoices      spatial members per octave of wavelength (1: octave bands)
%   ChannelMinTile   the position diagonal, part 1: per-CHANNEL rows are stored only at tiles
%                    at least this long (16 s); finer tiles are described by the sensor tree's
%                    nodes, and per-sensor detail below it is what phase 2 fetches. 0 keeps
%                    channels at every level. Ignored when the store has no sensor tree.
%   ChannelEnvelope  keep per-channel MIN and MAX at every level even below ChannelMinTile
%                    (true). Those two merge by extremum, so the pair at any level is a
%                    bound that CONTAINS every sample of its tile: that is a min/max mipmap,
%                    which is what draws an exact level-of-detail envelope of a trace at any
%                    zoom without touching raw samples. Measured on the reference data: 10.4 MB for the
%                    whole pyramid against a 128 MB store, because the diagonal exists for
%                    the band arrays, not for the all-pass moments. The other channel columns
%                    still start at ChannelMinTile.
%   Moments          "extended" (default) adds four more per-tile sums -- sum|x|, sum sqrt|x|,
%                    sum x^3, sum x^4 -- which merge exactly like the others and between them
%                    derive the whole shape family a time-domain feature extractor reports:
%                    skewness, kurtosis, shape factor, impulse factor, clearance factor.
%                    "basic" keeps only n, sumX, sumX2, absMax, min, max. Measured on the reference data:
%                    +33 MB on a 139 MB store, because the group rows carry them too.
%   Peaks            also store, for every band AT ITS OWN LEVEL, the spectral peak of each
%                    tile by FFT: frequency, amplitude and the band's power from that same
%                    spectrum (true; rheome.ingest.peaks). NOT mergeable -- and it does not need to
%                    be, because the diagonal gives each band exactly one home level. The
%                    bank's `energy` stays the exact power; this says where inside the band
%                    it sits, to a fraction of an FFT bin.
%   Couple           phase-amplitude coupling accumulators per tile: "auto" (every octave
%                    pair at least two octaves apart), false, or an [nP x 2] list of
%                    (slow, fast) band ids. Stored at the SLOW band's level, where the
%                    tiling gives 22.6 of its cycles in every band -- the same number, so
%                    the estimator's bias does not vary down the ladder. ⭐ THEY MERGE: the
%                    complex sum of A_fast*exp(i*phi_slow) and the sum of A_fast are sums,
%                    so the vector length and the preferred phase come back exactly at any
%                    coarser level (rheome.ingest.couple, rheome.select.coupling).
%   PhaseBins        also accumulate the modulogram: the amplitude sum per slow-phase bin
%                    (0 = no). Sums too, so Tort's modulation index merges as well; costs
%                    PhaseBins/3 times the coupling arrays, hence off by default.
%   SpaceDiagonal    the position diagonal, part 2: a wavelength band is stored only on nodes
%                    whose diameter is at least the band's shortest wavelength, plus the root
%                    (true). Needs calibrated wavelengths; otherwise every node carries every band.
%   MaxSupport       longest wavelet support carried, in seconds (Inf = whatever fits the
%                    record). The floor of the bank is support*fc / MaxSupport (6.46 s/Hz
%                    for Morse at 10 voices): 6.5 s -> 1 Hz, 13 s -> 0.5 Hz, 65 s -> 0.1 Hz.
%                    Cuts the transform (each scale is a full-rate pass) and bounds the
%                    page halo; the energy below the floor stays in sumX2 and the residual.
% OUTPUT:
%   cfg  .FrameFloor .Wavelet .VoicesPerOctave .FrequencyLimits .Precision .MaxBytes
%        .Paged .PageOverhead .MinPageLength .MaxSupport .Bank .Anchor .Oversample .Space
%        .SpaceVoices .ChannelMinTile .ChannelEnvelope .Moments .Peaks .Couple .PhaseBins
%        .SpaceDiagonal .hash
%
% ⚠ The page policy (Paged, PageOverhead, MinPageLength) is part of the hash: a paged
% store differs from a whole-record one inside the cone (design §3.12). MaxBytes is NOT:
% it is this machine's budget, and it changes only which job computes a band, never the
% band's values. Whether "auto" actually paged is recorded in meta.paged and meta.plan.
%
% See also: rheome.ingest.grid, rheome.ingest.bank, rheome.ingest.reduce
%
% Author: Diellor Basha, 2026

    arguments
        opts.FrameFloor      (1,1) double {mustBePositive} = 0.25
        opts.Wavelet         (1,1) string {mustBeMember(opts.Wavelet, ["morse","amor","bump"])} = "morse"
        opts.VoicesPerOctave       double = []
        opts.FrequencyLimits       double = []
        opts.Precision       (1,1) string {mustBeMember(opts.Precision, ["single","double"])} = "single"
        opts.MaxBytes        (1,1) double {mustBePositive} = 2.5e9
        opts.Paged           (1,1) string {mustBeMember(opts.Paged, ["auto","always","never"])} = "auto"
        opts.PageOverhead    (1,1) double {mustBePositive, mustBeLessThanOrEqual(opts.PageOverhead, 2)} = 0.25
        opts.MinPageLength   (1,1) double {mustBePositive} = 16
        opts.MaxSupport      (1,1) double {mustBePositive} = Inf
        opts.Bank            (1,1) string {mustBeMember(opts.Bank, ["frame","morse"])} = "frame"
        opts.Anchor          (1,1) double {mustBePositive} = 1
        opts.Oversample      (1,1) double {mustBeGreaterThanOrEqual(opts.Oversample, 1)} = 2
        opts.Space           (1,1) logical = true
        opts.SpaceVoices     (1,1) double {mustBeInteger, mustBePositive} = 1
        opts.ChannelMinTile  (1,1) double {mustBeNonnegative} = 16
        opts.ChannelEnvelope (1,1) logical = true
        opts.Moments         (1,1) string {mustBeMember(opts.Moments, ["basic","extended"])} = "extended"
        opts.Peaks           (1,1) logical = true
        opts.Couple                 = "auto"
        opts.PhaseBins       (1,1) double {mustBeInteger, mustBeNonnegative} = 0
        opts.SpaceDiagonal   (1,1) logical = true
    end

    if isempty(opts.VoicesPerOctave)
        % the frame bank's members are compactly supported in frequency, so their time
        % support grows with voices (measured: support*fc = 8, 15, 37 at 2, 4, 10 voices
        % against Morse's 2.9 under the same 99.9%-energy definition). Four keeps the alpha
        % band at 2 s tiles with Q ~ 6; Morse keeps the flow pipeline's ten.
        if opts.Bank == "frame", opts.VoicesPerOctave = 4; else, opts.VoicesPerOctave = 10; end
    elseif ~(isscalar(opts.VoicesPerOctave) && opts.VoicesPerOctave >= 1 && opts.VoicesPerOctave == round(opts.VoicesPerOctave))
        error('ingest:config:voices', 'VoicesPerOctave must be a positive integer.');
    end
    if ~isempty(opts.FrequencyLimits)
        fl = opts.FrequencyLimits(:)';
        if numel(fl) ~= 2 || ~(fl(1) > 0) || ~(fl(2) > fl(1))
            error('ingest:config:limits', 'FrequencyLimits must be [] or [fLo fHi] with 0 < fLo < fHi.');
        end
        opts.FrequencyLimits = fl;
    end

    % canonical order -- assignment order IS field order, and the hash depends on it
    cfg = struct();
    cfg.FrameFloor      = opts.FrameFloor;
    cfg.Wavelet         = char(opts.Wavelet);
    cfg.VoicesPerOctave = opts.VoicesPerOctave;
    cfg.FrequencyLimits = opts.FrequencyLimits;
    cfg.Precision       = char(opts.Precision);
    cfg.MaxBytes        = opts.MaxBytes;
    cfg.Paged           = char(opts.Paged);
    cfg.PageOverhead    = opts.PageOverhead;
    cfg.MinPageLength   = opts.MinPageLength;
    cfg.MaxSupport      = opts.MaxSupport;
    cfg.Bank            = char(opts.Bank);
    cfg.Anchor          = opts.Anchor;
    cfg.Oversample      = opts.Oversample;
    cfg.Space           = opts.Space;
    cfg.SpaceVoices     = opts.SpaceVoices;
    cfg.ChannelMinTile  = opts.ChannelMinTile;
    cfg.ChannelEnvelope = opts.ChannelEnvelope;
    cfg.Moments         = opts.Moments;
    cfg.Peaks           = opts.Peaks;
    cfg.Couple          = opts.Couple;
    cfg.PhaseBins       = opts.PhaseBins;
    cfg.SpaceDiagonal   = opts.SpaceDiagonal;
    h = rmfield(cfg, 'MaxBytes');                       % a machine budget, not a property of the values
    cfg.hash            = i_md5(jsonencode(h));
end

function h = i_md5(s)
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(s));
    d = typecast(md.digest(), 'uint8');
    h = lower(reshape(dec2hex(d, 2)', 1, []));
    h = h(1:8);
end
% Author: Diellor Basha, 2026
