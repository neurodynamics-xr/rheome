function disp(obj)
% DISP  One-screen summary: what the bank is, and which members to distrust.
%
% Author: Diellor Basha, 2026

    if ~isscalar(obj)
        fprintf('  %s array with size %s\n', class(obj), mat2str(size(obj)));
        return;
    end

    b  = framebounds(obj);
    nu = nnz(~obj.Usable);
    kk = centerWavenumbers(obj);
    fprintf('  graphfilterbank: %s, %d members\n', obj.Wavelet, obj.NumMembers);
    fprintf('    lambda_max %.4g  |  k range %.3g - %.3g rad/m\n', ...
        obj.Lmax_, min(kk), max(kk));
    fprintf('    sizes      %.4g - %.4g  |  %.1f voices/octave\n', ...
        obj.SizeLimits(1), obj.SizeLimits(2), obj.AchievedVoicesPerOctave);
    fprintf('    frame      A = %.4g  B = %.4g  B/A = %.4g%s\n', ...
        b.A, b.B, b.Tightness, i_tern(isframetight(obj), '  (tight)', ''));
    if nu > 0
        fprintf('    %d unusable member(s): >5%% of their mass lies beyond lambda_max,\n', nu);
        fprintf('       so they under-respond. Usable floor is sigma >= %.4g.\n', obj.FinestScale);
    else
        fprintf('    all members usable (floor sigma >= %.4g)\n', obj.FinestScale);
    end
    if isempty(obj.Transform)
        fprintf('    no Transform attached: design queries only.\n');
        fprintf('       attach rheome.graphtransform.eigen(Phi,B,Lambda) to filter fields.\n');
    end
end

function s = i_tern(c, a, b)
    if c, s = a; else, s = b; end
end

% Author: Diellor Basha, 2026
