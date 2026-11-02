function r = evenround(x)
% EVENROUND  Round half to even, as Python's round and numpy's rint do (MATLAB's round goes away from zero).
    r = round(x);
    tie = abs(x - fix(x)) == 0.5;
    r(tie) = 2 * round(x(tie) / 2);
end
