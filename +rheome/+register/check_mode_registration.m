function C = check_mode_registration(T, R, subject, template, reg)
% REGISTER.CHECK_MODE_REGISTRATION  A registration in the modes measured on smooth fields: the real spherical harmonics l = 1..3.
%
%   C = rheome.register.check_mode_registration(T, R, subject, template, reg)        % [T, R] = rheome.register.mode_registration(...)
%
% C  table, one row per (l, m):
%    errReg    the registration: R's pull of the subject's field, analysed on the template, against the field analysed there
%    errTrunc  T against that pull: what the subject's band loses
%    errT      T's result against the field analysed on the template;  corrT  the synthesised result against the field
%              (area-weighted correlation on the template's glued vertices)
%    gate      errReg <= 0.05;  target  errT < 0.15 and corrT > 0.95
% Each error is the norm of the difference over the norm of the truth (coefficients), so the harmonics' normalisation does not
% enter. As nxr.register.modes.check_mode_registration (err_reg, err_trunc, err_T, corr_T).
%
% See also: rheome.register.mode_registration
%
% Author: Diellor Basha, 2026
arguments, T (:,:) double, R (:,:), subject (1,1) struct, template (1,1) struct, reg (1,1) struct, end
ms = subject.kept; mt = template.kept;
[Ys, lm] = i_harmonics(reg.subject_sphere(subject.keep(1:ms), :));
Ys = [Ys; zeros(size(subject.positions, 1) - ms, size(Ys, 2))];  % the subject's bridge: not cortex, no data
[cut, B] = bridge_rows(template, reg.template_positions);
Yt = [i_harmonics(reg.template_sphere(template.keep(1:mt), :)); B * i_harmonics(reg.template_sphere(cut, :))];
cs = fold(subject, Ys); ct = fold(template, Yt); creg = fold(template, R * Ys);
Tc = T * cs; rec = double(template.P) * Tc;
nrm = @(x, y) vecnorm(x - y) ./ vecnorm(y);
m = template.mass(:);
x = rec - (m' * rec) / sum(m); y = Yt - (m' * Yt) / sum(m);
corrT = (sum(m .* x .* y) ./ sqrt(sum(m .* x .^ 2) .* sum(m .* y .^ 2)))';
errReg = nrm(creg, ct)'; errTrunc = nrm(Tc, creg)'; errT = nrm(Tc, ct)';
C = table(lm(:, 1), lm(:, 2), errReg, errTrunc, errT, corrT, errReg <= 0.05, errT < 0.15 & corrT > 0.95, ...
    'VariableNames', {'l', 'm', 'errReg', 'errTrunc', 'errT', 'corrT', 'gate', 'target'});
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
