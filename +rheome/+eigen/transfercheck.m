function C = transfercheck(T, R, subj, tmpl)
% EIGEN.TRANSFERCHECK  A transfer measured on smooth fields: the real spherical harmonics l = 1..3 of the sphere positions.
%
%   C = rheome.eigen.transfercheck(T, R, subj, tmpl)        % [T, R] = rheome.eigen.transfer(subj, tmpl)
%
% C  table, one row per (l, m):
%    errReg    the registration: R's pull of the subject's field, analysed on the template, against the field analysed there
%    errTrunc  T against that pull: what the subject's band loses
%    errT      T's result against the field analysed on the template;  corrT  the synthesised result against the field
%              (area-weighted correlation on the template's glued vertices)
%    gate      errReg <= 0.05;  target  errT < 0.15 and corrT > 0.95
% Each error is the norm of the difference over the norm of the truth (coefficients), so the harmonics' normalisation does not
% enter. Port of nxr-convert eigen/register.py check_transfer.
%
% See also: rheome.eigen.transfer
%
% Author: Diellor Basha, 2026 (port of nxr-convert check_transfer)
arguments, T (:,:) double, R (:,:), subj (1,1) struct, tmpl (1,1) struct, end
ms = nnz(subj.keep > 0); mt = nnz(tmpl.keep > 0);
[Ys, lm] = i_harmonics(subj.Sphere(subj.keep(1:ms), :));
Ys = [Ys; zeros(size(subj.X, 1) - ms, size(Ys, 2))];             % the subject's bridge: not cortex, no data
cut = setdiff((1:size(tmpl.Pin, 1))', tmpl.keep(1:mt));
[nn, d] = knnsearch(tmpl.Pin(cut, :), tmpl.X(mt+1:end, :), 'K', 3);
w = 1 ./ max(d, 1e-9); w = w ./ sum(w, 2);
Yc = i_harmonics(tmpl.Sphere(cut, :));
Yb = zeros(size(tmpl.X, 1) - mt, size(Yc, 2)); for j = 1:3, Yb = Yb + w(:, j) .* Yc(nn(:, j), :); end
Yt = [i_harmonics(tmpl.Sphere(tmpl.keep(1:mt), :)); Yb];
cs = i_fold(subj.pkg, Ys); ct = i_fold(tmpl.pkg, Yt); creg = i_fold(tmpl.pkg, R * Ys);
Tc = T * cs; rec = double(tmpl.pkg.P) * Tc;
nrm = @(x, y) vecnorm(x - y) ./ vecnorm(y);
m = tmpl.pkg.mass(:);
x = rec - (m' * rec) / sum(m); y = Yt - (m' * Yt) / sum(m);
corrT = (sum(m .* x .* y) ./ sqrt(sum(m .* x .^ 2) .* sum(m .* y .^ 2)))';
errReg = nrm(creg, ct)'; errTrunc = nrm(Tc, creg)'; errT = nrm(Tc, ct)';
C = table(lm(:, 1), lm(:, 2), errReg, errTrunc, errT, corrT, errReg <= 0.05, errT < 0.15 & corrT > 0.95, ...
    'VariableNames', {'l', 'm', 'errReg', 'errTrunc', 'errT', 'corrT', 'gate', 'target'});
end

function E = i_fold(pkg, X)
P = double(pkg.P); co = pkg.kind(:) == 2;
E = P' * (pkg.mass(:) .* X); E(~co, :) = 2 * E(~co, :);
if any(co), G = P(:, co)' * (pkg.mass(:) .* P(:, co)); Rc = chol((G + G') / 2); E(co, :) = Rc \ (Rc' \ E(co, :)); end
end

function [Y, lm] = i_harmonics(X)
% real spherical harmonics l = 1..3 of unit positions: cos(m phi) for m >= 0, sin(|m| phi) for m < 0, times P_l^|m|(cos theta)
X = X ./ vecnorm(X, 2, 2);
th = acos(max(-1, min(1, X(:, 3)))); ph = atan2(X(:, 2), X(:, 1));
Y = []; lm = [];
for l = 1:3
    Pl = legendre(l, cos(th), 'norm')';                          % [n x (l+1)], orders 0..l
    for m = -l:l
        if m < 0, y = Pl(:, -m + 1) .* sin(-m * ph); else, y = Pl(:, m + 1) .* cos(m * ph); end
        Y = [Y y]; lm = [lm; l m]; %#ok<AGROW>
    end
end
end
