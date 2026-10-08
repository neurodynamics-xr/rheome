function A = level(db, L, name)
% SELECT.LEVEL  One level array from the store, cached, with its bytes counted.
%
%   A = rheome.select.level(db, L, name)     name in n, sumX, sumX2, absMax, min, max, energy, envMax, nCoi
%
% Author: Diellor Basha, 2026

    key = sprintf('L%02d_%s', L, name);
    if isKey(db.cache, key)
        A = db.cache(key);
        return
    end
    A = db.m.(key);
    db.cache(key) = A;
    w = whos('A');
    db.cost('bytes') = db.cost('bytes') + w.bytes;
    db.cost('reads') = db.cost('reads') + 1;
end
% Author: Diellor Basha, 2026
