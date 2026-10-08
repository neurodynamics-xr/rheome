function T = ladder(db, opts)
% SELECT.LADDER  The constant-Q ladder of a store: what each band costs in time and in rate.
%
%   T = rheome.select.ladder(db)
%   T = rheome.select.ladder(db, Scope="channel")     % adds the floor a single sensor imposes
%
% One row per band, answering the question a level-of-detail view keeps asking: *at this
% band, how long is a tile, how many cycles is that, and at what rate does the band arrive?*
%
% ⭐ THE TWO INVARIANTS ARE THE POINT. A band's support is inversely proportional to its
% centre frequency (constant Q), and its tile is the first dyadic level at least that long,
% so **cycles per tile is the same in every octave** — 22.6 on the reference bank. Each member
% is evaluated at a rate proportional to its centre frequency (0.7*fc here), so **samples
% per tile is the same in every octave** too — 19. An octave of the spectrum therefore costs
% the same to describe wherever it sits, while the record's own samples per tile double with
% every level: 150 at 64-128 Hz, 614 400 at the bottom band. That ratio is the whole reason
% for a multirate constant-Q store rather than a spectrogram.
%
% COLUMNS
%   band_id, kind, f_lo, f_hi, f_center, octaves, q
%   support_s        the member's time support (99.9 % of energy), seconds
%   natural_level    the first level whose tile is at least the support
%   tile_s           that level's tile length
%   cycles_per_tile  tile_s * f_center
%   tiles_in_record  how many such tiles the record holds
%   rate_hz          the rate the band's members are evaluated at (NaN for a Morse store)
%   samples_per_tile tile_s * rate_hz, the band's own cost per tile
%   samples_full     tile_s * fs, what the record costs for the same tile
%   floor_level, floor_tile_s   with Scope="channel", the coarser of the natural level and
%                    grid.channelLevel: what a SINGLE SENSOR can actually be read at
%
% See also: rheome.select.derive, rheome.recordingbrowser, rheome.ingest.bank, rheome.timefilterbank/rates
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        opts.Scope (1,1) string {mustBeMember(opts.Scope, ["channel","group","none"])} = "none"
    end
    b = db.bands;  g = db.grid;  m = db.meta;  nB = height(b);
    L = b.naturalLevel;
    tile = g.tExtent(L + 1)';
    r = i_rates(db);

    T = table(b.j, string(b.kind), b.fLo, b.fHi, b.fCenter, b.fExtent, b.fCenter ./ (b.fHi - b.fLo), ...
              b.tSupport, L, tile, tile .* b.fCenter, g.K(L + 1)', ...
              r, tile .* r, tile * m.fs, ...
              'VariableNames', {'band_id','kind','f_lo','f_hi','f_center','octaves','q', ...
                                'support_s','natural_level','tile_s','cycles_per_tile', ...
                                'tiles_in_record','rate_hz','samples_per_tile','samples_full'});
    if opts.Scope == "channel"
        fl = max(L, db.channelLevel);
        T.floor_level = fl;
        T.floor_tile_s = g.tExtent(fl + 1)';
    end
end

% The rate comes from the bank's design, not from transforming anything: rebuilding a
% timefilterbank for the record's length costs 60 ms and is exact (@timefilterbank/rates).
% A Morse store's members are not evaluated per band that way, so it reports NaN rather
% than a number that would mean something else.
function r = i_rates(db)
    m = db.meta;  cfg = m.cfg;  nB = height(db.bands);
    r = nan(nB, 1);
    if ~isfield(cfg, 'Bank') || ~strcmpi(char(cfg.Bank), 'frame'), return; end
    fb = rheome.timefilterbank(m.nT, 'SamplingFrequency', m.fs, ...
                        'VoicesPerOctave', cfg.VoicesPerOctave, 'Anchor', cfg.Anchor, ...
                        'FrequencyLimits', cfg.FrequencyLimits, 'Oversample', cfg.Oversample);
    mr = rates(fb);
    for j = 1:nB
        mem = db.bands.scales{j};
        if ~isempty(mem), r(j) = max(mr(mem)); end       % the band's fastest member
    end
end
% Author: Diellor Basha, 2026
