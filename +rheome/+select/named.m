function T = named(db, opts)
% SELECT.NAMED  The physiological bands, mapped onto the store's octave bands.
%
%   T = rheome.select.named(db)
%   T = rheome.select.named(db, Bands=struct('alpha', [8 13], 'beta', [13 30]))
%
% ⭐ WHAT THIS IS AND IS NOT. The store's bands are OCTAVES on an absolute grid anchored at
% 1 Hz, because that is what makes them constant-Q, mergeable and comparable across
% subjects and sampling rates. Delta, theta, alpha, beta and gamma are conventions with
% ragged edges that no dyadic grid can hit. This table says, for each name, which stored
% bands cover it, how well, and at which level they live -- so a query can say "alpha" and
% get an answer whose approximation is on the record instead of in someone's head.
%
% ⚠ READ `covered` BEFORE QUOTING A NUMBER. alpha 8-13 Hz is covered by the 8-16 Hz octave,
% which is 1.6x too wide: 62 % of the octave's width is alpha. beta 13-30 Hz straddles two
% octaves and is covered to 89 %. Where the fit is poor, the honest answer is the octave's
% own edges, which is why they are in the table too.
%
% COLUMNS
%   name, f_lo, f_hi       the convention
%   band_ids               the stored octave bands that overlap it
%   band_lo, band_hi       what those octaves actually span
%   covered                the fraction of the octave span that is inside the convention
%   natural_level, tile_s  where those bands live, and their tile length
%   peak_level             the level whose stored spectral peak covers this name
%
% See also: rheome.select.ladder, rheome.select.derive, rheome.ingest.peaks
%
% Author: Diellor Basha, 2026

    arguments
        db (1,1) struct
        opts.Bands = struct('delta', [1 4], 'theta', [4 8], 'alpha', [8 13], ...
                            'beta', [13 30], 'gamma', [30 80])
    end
    b = db.bands;
    oct = strcmp(string(b.kind), "octave");
    nm = fieldnames(opts.Bands);
    name = strings(0,1); flo = []; fhi = []; ids = {}; blo = []; bhi = [];
    cov = []; lev = {}; tile = {}; pk = [];
    for i = 1:numel(nm)
        r = opts.Bands.(nm{i});
        k = find(oct & b.fHi > r(1) & b.fLo < r(2));
        if isempty(k), continue; end
        name(end+1,1) = string(nm{i});                     %#ok<AGROW>
        flo(end+1,1) = r(1);  fhi(end+1,1) = r(2);         %#ok<AGROW>
        ids{end+1,1} = b.j(k)';                            %#ok<AGROW>
        blo(end+1,1) = min(b.fLo(k));  bhi(end+1,1) = max(b.fHi(k));   %#ok<AGROW>
        ov = max(0, min(bhi(end), r(2)) - max(blo(end), r(1)));
        cov(end+1,1) = ov / (bhi(end) - blo(end));         %#ok<AGROW>
        lev{end+1,1} = b.naturalLevel(k)';                 %#ok<AGROW>
        tile{end+1,1} = db.grid.tExtent(b.naturalLevel(k) + 1);        %#ok<AGROW>
        pk(end+1,1) = min(b.naturalLevel(k));              %#ok<AGROW>
    end
    T = table(name, flo, fhi, ids, blo, bhi, cov, lev, tile, pk, ...
              'VariableNames', {'name','f_lo','f_hi','band_ids','band_lo','band_hi', ...
                                'covered','natural_level','tile_s','peak_level'});
end
% Author: Diellor Basha, 2026
