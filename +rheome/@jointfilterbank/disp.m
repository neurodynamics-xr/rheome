function disp(obj)
% DISP  One-screen summary: shape, coverage, and what materialising would cost.
% Author: Diellor Basha, 2026

    if ~isscalar(obj)
        fprintf('  %s array with size %s\n', class(obj), mat2str(size(obj)));
        return;
    end
    K = numel(obj.Lambda);  nO = obj.NumOmega;  Nf = obj.NumMembers;
    b = framebounds(obj);
    bytes = K * nO * Nf * 8;

    fprintf('  jointfilterbank: %d members = %d graph x %d time x %d joint%s\n', ...
        Nf, obj.NumGraph, obj.NumTime, obj.NumJoint, ...
        i_tern(obj.Separable, '  (separable)', '  (non-separable)'));
    fprintf('    grid       %d modes x %d bins  |  %.4g - %.4g Hz  |  half ''%s'', boundary ''%s''\n', ...
        K, nO, min(obj.Frequencies), max(obj.Frequencies), obj.Half, obj.Boundary);
    fprintf('    frame      A = %.4g  B = %.4g  B/A = %.4g  |  %.0f%% uncovered\n', ...
        b.A, b.B, b.Tightness, 100*b.Uncovered);
    fprintf('    materialised the bank would be %s -- stored as factors instead\n', jfb_bytes(bytes));
    if isempty(obj.Transform)
        fprintf('    no Transform attached: joint-domain queries only.\n');
        fprintf('       attach rheome.graphtransform.eigen(Phi,B,Lambda) for impulse responses.\n');
    end
end

function s = i_tern(c, a, b)
    if c, s = a; else, s = b; end
end

% Author: Diellor Basha, 2026
