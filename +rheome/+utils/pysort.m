function a = pysort(k)
% UTILS.PYSORT  CPython 3.12's list.sort (TimSort, powersort merge policy, galloping) on 1:numel(k) with
% the comparator x < y  <=>  k(x) - k(y) < 0 (NaN reads as "not less"), step for step.
%
% Why: nsp's atlas ladder sorts a disconnected tile of >= 64 vertices with Python's sorted() and a
% comparator that is inconsistent where keys are NaN; the result then depends on the exact steps.
% A port of Objects/listobject.c (count_run, binarysort, powerloop, merge_at, gallop_left/right,
% merge_lo/hi, merge_force_collapse). Pinned against Python 3.12 by tests/tSurfaceLadder.m.
%
% Author: Diellor Basha, 2026

    k = k(:);  n = numel(k);  a = (1:n)';
    if n < 2, return; end
    lt = @(x, y) k(x) - k(y) < 0;
    minGallop = 7;  MIN_GALLOP = 7;
    runBase = zeros(0, 1);  runLen = runBase;  runPow = runBase;     % pending runs (0-based bases)
    m = n;  r = 0;
    while m >= 64, r = bitor(r, bitand(m, 1));  m = bitshift(m, -1); end
    minrun = m + r;

    lo = 0;  rem = n;
    while rem > 0
        % count_run
        if rem == 1
            nr = 1;  desc = false;
        else
            nr = 2;  desc = lt(a(lo + 2), a(lo + 1));
            for j = lo + 3 : lo + rem
                if lt(a(j), a(j - 1)) ~= desc, break; end
                nr = nr + 1;
            end
        end
        if desc, a(lo + 1 : lo + nr) = a(lo + nr : -1 : lo + 1); end
        if nr < minrun
            force = min(rem, minrun);
            binarysort(lo, lo + force, lo + nr);
            nr = force;
        end
        % found_new_run
        if ~isempty(runLen)
            pw = powerloop(runBase(end), runLen(end), nr, n);
            while numel(runLen) > 1 && runPow(end - 1) > pw
                mergeAt(numel(runLen) - 1);
            end
            runPow(end) = pw;
        end
        runBase(end + 1, 1) = lo;  runLen(end + 1, 1) = nr;  runPow(end + 1, 1) = 0;
        lo = lo + nr;  rem = rem - nr;
    end
    while numel(runLen) > 1                                         % merge_force_collapse
        i = numel(runLen) - 1;
        if i > 1 && runLen(i - 1) < runLen(i + 1), i = i - 1; end
        mergeAt(i);
    end

    % sorts a(lo+1 : hi), a(lo+1 : start) already sorted (0-based half-open bounds)
    function binarysort(lo, hi, start)
        if lo == start, start = start + 1; end
        for s = start : hi - 1
            pv = a(s + 1);  l = lo;  rr = s;
            while l < rr
                p = l + bitshift(rr - l, -1);
                if lt(pv, a(p + 1)), rr = p; else, l = p + 1; end
            end
            a(l + 2 : s + 1) = a(l + 1 : s);  a(l + 1) = pv;
        end
    end

    % merge pending runs i and i+1 (1-based stack index)
    function mergeAt(i)
        sa = runBase(i);  na = runLen(i);  sb = runBase(i + 1);  nb = runLen(i + 1);
        runLen(i) = na + nb;
        runBase(i + 1) = [];  runLen(i + 1) = [];  runPow(i + 1) = [];   % p[i+1] = p[i+2], --n
        kk = gallopRight(a(sb + 1), a(sa + 1 : sa + na), 0);
        sa = sa + kk;  na = na - kk;
        if na == 0, return; end
        nb = gallopLeft(a(sa + na), a(sb + 1 : sb + nb), nb - 1);
        if nb <= 0, return; end
        if na <= nb, mergeLo(sa, na, sb, nb); else, mergeHi(sa, na, sb, nb); end
    end

    function mergeLo(sa, na, sb, nb)
        A = a(sa + 1 : sa + na);  ia = 1;  dest = sa;                % A = temp copy of run a
        dest = dest + 1;  a(dest) = a(sb + 1);  sb = sb + 1;  nb = nb - 1;
        if nb == 0, finishLo(); return; end
        if na == 1, copyB(); return; end
        mg = minGallop;
        while true
            ac = 0;  bc = 0;
            while true
                if lt(a(sb + 1), A(ia))
                    dest = dest + 1;  a(dest) = a(sb + 1);  sb = sb + 1;  bc = bc + 1;  ac = 0;  nb = nb - 1;
                    if nb == 0, finishLo(); return; end
                    if bc >= mg, break; end
                else
                    dest = dest + 1;  a(dest) = A(ia);  ia = ia + 1;  ac = ac + 1;  bc = 0;  na = na - 1;
                    if na == 1, copyB(); return; end
                    if ac >= mg, break; end
                end
            end
            mg = mg + 1;
            while true
                mg = mg - (mg > 1);  minGallop = mg;
                kk = gallopRight(a(sb + 1), A(ia : ia + na - 1), 0);  ac = kk;
                if kk
                    a(dest + 1 : dest + kk) = A(ia : ia + kk - 1);  dest = dest + kk;  ia = ia + kk;  na = na - kk;
                    if na == 1, copyB(); return; end
                    if na == 0, finishLo(); return; end
                end
                dest = dest + 1;  a(dest) = a(sb + 1);  sb = sb + 1;  nb = nb - 1;
                if nb == 0, finishLo(); return; end
                kk = gallopLeft(A(ia), a(sb + 1 : sb + nb), 0);  bc = kk;
                if kk
                    a(dest + 1 : dest + kk) = a(sb + 1 : sb + kk);  dest = dest + kk;  sb = sb + kk;  nb = nb - kk;
                    if nb == 0, finishLo(); return; end
                end
                dest = dest + 1;  a(dest) = A(ia);  ia = ia + 1;  na = na - 1;
                if na == 1, copyB(); return; end
                if ~(ac >= MIN_GALLOP || bc >= MIN_GALLOP), break; end
            end
            mg = mg + 1;  minGallop = mg;
        end
        function finishLo()
            if na, a(dest + 1 : dest + na) = A(ia : ia + na - 1); end
        end
        function copyB()
            a(dest + 1 : dest + nb) = a(sb + 1 : sb + nb);  a(dest + nb + 1) = A(ia);
        end
    end

    function mergeHi(sa, na, sb, nb)
        B = a(sb + 1 : sb + nb);  ib = nb;                           % B = temp copy of run b
        dest = sb + nb;  pa = sa + na;                               % 1-based: current a is a(pa), dest a(dest)
        a(dest) = a(pa);  dest = dest - 1;  pa = pa - 1;  na = na - 1;
        if na == 0, finishHi(); return; end
        if nb == 1, copyA(); return; end
        mg = minGallop;
        while true
            ac = 0;  bc = 0;
            while true
                if lt(B(ib), a(pa))
                    a(dest) = a(pa);  dest = dest - 1;  pa = pa - 1;  ac = ac + 1;  bc = 0;  na = na - 1;
                    if na == 0, finishHi(); return; end
                    if ac >= mg, break; end
                else
                    a(dest) = B(ib);  dest = dest - 1;  ib = ib - 1;  bc = bc + 1;  ac = 0;  nb = nb - 1;
                    if nb == 1, copyA(); return; end
                    if bc >= mg, break; end
                end
            end
            mg = mg + 1;
            while true
                mg = mg - (mg > 1);  minGallop = mg;
                kk = gallopRight(B(ib), a(sa + 1 : sa + na), na - 1);
                kk = na - kk;  ac = kk;
                if kk
                    a(dest - kk + 1 : dest) = a(pa - kk + 1 : pa);  dest = dest - kk;  pa = pa - kk;  na = na - kk;
                    if na == 0, finishHi(); return; end
                end
                a(dest) = B(ib);  dest = dest - 1;  ib = ib - 1;  nb = nb - 1;
                if nb == 1, copyA(); return; end
                kk = gallopLeft(a(pa), B(1 : nb), nb - 1);
                kk = nb - kk;  bc = kk;
                if kk
                    a(dest - kk + 1 : dest) = B(ib - kk + 1 : ib);  dest = dest - kk;  ib = ib - kk;  nb = nb - kk;
                    if nb == 1, copyA(); return; end
                    if nb == 0, finishHi(); return; end
                end
                a(dest) = a(pa);  dest = dest - 1;  pa = pa - 1;  na = na - 1;
                if na == 0, finishHi(); return; end
                if ~(ac >= MIN_GALLOP || bc >= MIN_GALLOP), break; end
            end
            mg = mg + 1;  minGallop = mg;
        end
        function finishHi()
            if nb, a(dest - nb + 1 : dest) = B(1 : nb); end
        end
        function copyA()
            a(dest - na + 1 : dest) = a(pa - na + 1 : pa);  dest = dest - na;  a(dest) = B(ib);
        end
    end

    % k such that v(k) < key <= v(k+1) (0-based k into the run v), starting at hint
    function kk = gallopLeft(key, v, hint)
        nn = numel(v);  last = 0;  ofs = 1;
        if lt(v(hint + 1), key)
            maxofs = nn - hint;
            while ofs < maxofs
                if lt(v(hint + ofs + 1), key), last = ofs;  ofs = ofs * 2 + 1; else, break; end
            end
            ofs = min(ofs, maxofs);  last = last + hint;  ofs = ofs + hint;
        else
            maxofs = hint + 1;
            while ofs < maxofs
                if lt(v(hint - ofs + 1), key), break; end
                last = ofs;  ofs = ofs * 2 + 1;
            end
            ofs = min(ofs, maxofs);  t = last;  last = hint - ofs;  ofs = hint - t;
        end
        last = last + 1;
        while last < ofs
            mm = last + bitshift(ofs - last, -1);
            if lt(v(mm + 1), key), last = mm + 1; else, ofs = mm; end
        end
        kk = ofs;
    end

    % k such that v(k) <= key < v(k+1) (0-based), starting at hint
    function kk = gallopRight(key, v, hint)
        nn = numel(v);  last = 0;  ofs = 1;
        if lt(key, v(hint + 1))
            maxofs = hint + 1;
            while ofs < maxofs
                if lt(key, v(hint - ofs + 1)), last = ofs;  ofs = ofs * 2 + 1; else, break; end
            end
            ofs = min(ofs, maxofs);  t = last;  last = hint - ofs;  ofs = hint - t;
        else
            maxofs = nn - hint;
            while ofs < maxofs
                if lt(key, v(hint + ofs + 1)), break; end
                last = ofs;  ofs = ofs * 2 + 1;
            end
            ofs = min(ofs, maxofs);  last = last + hint;  ofs = ofs + hint;
        end
        last = last + 1;
        while last < ofs
            mm = last + bitshift(ofs - last, -1);
            if lt(key, v(mm + 1)), ofs = mm; else, last = mm + 1; end
        end
        kk = ofs;
    end
end

function r = powerloop(s1, n1, n2, n)
    r = 0;  a = 2 * s1 + n1;  b = a + n1 + n2;
    while true
        r = r + 1;
        if a >= n
            a = a - n;  b = b - n;
        elseif b >= n
            break
        end
        a = a * 2;  b = b * 2;
    end
end
