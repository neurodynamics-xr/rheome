function R = totable(file, opts)
% INGEST.TOTABLE  A store, or part of it, as long-form rows for a database or groupsummary.
%
%   R = rheome.ingest.totable(file)                              % top level, every stat
%   R = rheome.ingest.totable(file, Levels=[3 4], Channels=1:10, Bands=[], Stats=["energy","envMax"])
%
% One row per (channel, tile, band, stat). All-pass statistics carry band = 0 and NaN
% frequency extents. A band appears at a level only from its natural level up (the
% store's diagonal, grid.bandsAt); rows are keyed by the band's j, never by its slot. The row is the record the production database will hold, and the
% tile is described by centre + extent on both axes exactly as a query needs it
% (design §2, §3.10). Reads level variables through matfile, so only the requested levels
% touch disk.
%
% ⚠ Level 0 of a 600 s record at 270 channels and 15 bands is ~10 M rows per stat.
% Levels defaults to the TOP level for that reason; ask for finer levels explicitly.
%
% INPUTS:
%   file      store written by rheome.ingest.build
%   Levels    levels to emit ([] = top level only)
%   Channels  indices into the store's channel axis ([] = all)
%   Bands     band indices j ([] = all); all-pass stats are always included as band 0
%   Stats     subset of ["n","sumX","sumX2","absMax","min","max","energy","envMax","nCoi"]
% OUTPUT:
%   R  table: recording, channel, level, tCenter, tExtent, band, fCenter, fExtent, stat, value
%
% See also: rheome.ingest.build, rheome.flowfeatures/totable
%
% Author: Diellor Basha, 2026

    arguments
        file (1,:) char
        opts.Levels   double = []
        opts.Channels double = []
        opts.Bands    double = []
        opts.Stats    string = ["n","sumX","sumX2","absMax","min","max","energy","envMax","nCoi"]
    end

    m = matfile(file);
    meta = m.meta;  bands = m.bands;  g = m.grid;
    levels = opts.Levels;   if isempty(levels), levels = g.Lmax; end
    chans  = opts.Channels; if isempty(chans),  chans  = 1:meta.C; end
    bsel   = opts.Bands;    if isempty(bsel),   bsel   = (1:height(bands))'; end
    bsel = bsel(:);
    perRecord = ["n", "nCoi"];
    perBand   = ["energy", "envMax", "nCoi"];

    parts = {};
    for L = levels(:)'
        K = g.K(L+1);  tc = g.tCenter{L+1}(:);  te = g.tExtent(L+1);
        present = g.bandsAt{L+1};                       % bands stored at this level (diagonal)
        for st = opts.Stats
            v = m.(sprintf('L%02d_%s', L, st));
            if any(st == perBand)
                for b = bsel'
                    k = find(present == b, 1);
                    if isempty(k), continue; end        % below its natural level: not stored
                    if any(st == perRecord)
                        val = double(v(:, k));  ch = zeros(K, 1);
                    else
                        val = reshape(double(v(:, chans, k)), [], 1);
                        ch  = repelem(chans(:), K);
                    end
                    n = numel(val);
                    parts{end+1} = i_rows(meta.name, ch, L, repmat(tc, n/K, 1), te, b, ...
                                          bands.fCenter(b), bands.fExtent(b), st, val); %#ok<AGROW>
                end
            else
                if any(st == perRecord)
                    val = double(v(:, 1));  ch = zeros(K, 1);
                else
                    val = reshape(double(v(:, chans)), [], 1);
                    ch  = repelem(chans(:), K);
                end
                n = numel(val);
                parts{end+1} = i_rows(meta.name, ch, L, repmat(tc, n/K, 1), te, 0, NaN, NaN, st, val); %#ok<AGROW>
            end
        end
    end
    R = vertcat(parts{:});
end

function T = i_rows(name, ch, L, tc, te, b, fc, fe, st, val)
    n = numel(val);
    T = table(repmat(string(name), n, 1), ch(:), repmat(L, n, 1), tc(:), repmat(te, n, 1), ...
              repmat(b, n, 1), repmat(fc, n, 1), repmat(fe, n, 1), repmat(st, n, 1), val(:), ...
              'VariableNames', {'recording','channel','level','tCenter','tExtent', ...
                                'band','fCenter','fExtent','stat','value'});
end
% Author: Diellor Basha, 2026
