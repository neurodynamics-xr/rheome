function out = sensor_topology(outDir)
% DEMOS.SENSOR_TOPOLOGY  Divergence and curl of the phase gradient -- and what they classify.
%
%   rheome.demos.sensor_topology()
%   out = rheome.demos.sensor_topology('_figures')
%
%   1 field     the wavevector k = grad(phi) for four generators
%   2 div       divergence of k: where wavefronts are born
%   3 curl      curl of k, and the integer winding that is its honest form
%   4 classify  the two invariants, and the four families landing in four quadrants
%
% ⚠ A SCALAR SIGNAL HAS NO DIVERGENCE AND NO CURL. Those operators act on VECTOR fields, and
% one number per sensor is not one. The vector field that IS available is the local
% wavevector k = grad(phi) from the analytic phase, and the divergence and curl OF THAT are
% what this demo measures. Saying "the curl of the signal" is a category error; saying "the
% curl of its phase gradient" is a measurement.
%
% ⭐ CURL OF A GRADIENT IS ZERO -- WHICH IS EXACTLY WHY IT DETECTS ROTORS. curl(grad(phi))
% vanishes wherever phi is smooth and single-valued, so a curl map is ~0 everywhere EXCEPT
% where the phase is multivalued: a singularity. By Stokes the circulation of k around a loop
% is 2*pi times the winding number, so rheome.detect.phasesingularity's integer IS this curl in its
% integrated form -- and an integer cannot drift, where a density can.
%
% ⭐ THE TWO INVARIANTS SEPARATE THE FAMILIES, and this is the answer to "can it characterise
% a spiral":
%
%       generator     winding   divergence      because phi is
%       travelling      0          0            k.x        -- k constant
%       target          0        > 0 at focus   k r        -- k radial, fronts born at a point
%       rotating       +-1         0            m theta    -- k azimuthal, closed fronts
%       spiral         +-1       > 0            m theta + k r
%
% A spiral is precisely the pattern carrying BOTH. Winding alone cannot separate it from a
% rotor; divergence alone cannot separate it from a target. Together they do.
%
% ⚠ RUN THIS ON A FLAT ARRAY. Divergence of a tangent field on a CURVED array picks up a
% mean-curvature term (rheome.differential.divergence uses the full 3-vector deliberately), so a
% helmet reports divergence that is partly its own shape. A planar lattice has none.
%
% ⚠ |k| DIVERGES AT A CORE, so div and curl are ill-conditioned in the few faces around one.
% The maps are reported away from the core, and the winding -- which is exact there -- is what
% the core itself is measured with.
%
% INPUTS:
%   outDir  '' for screen, or a folder for PNGs
% OUTPUT:
%   out  .ok .checks .stats (one row per generator: winding, peak divergence)
%
% See also: rheome.demos.sensor_phase, rheome.flow.phasegradient, rheome.detect.phasesingularity,
%           rheome.differential.divergence, rheome.differential.curl
%
% Author: Diellor Basha, 2026

    if nargin < 1, outDir = ''; end
    doExport = ~isempty(outDir);
    chk = struct('name', {}, 'measured', {}, 'expected', {}, 'tol', {}, 'unit', {}, 'pass', {});
    fs = 200;  nT = 256;  tv = (0:nT-1)/fs;  f0 = 10;  jm = 128;

    % A flat, well-resolved lattice: divergence here is the field's, not the array's shape.
    arr = rheome.sensors.grid('Size', [24 24], 'Pitch', 10e-3);
    G   = rheome.sensors.graph(arr);
    S   = struct('Vertices', G.Vertices, 'Faces', G.Faces, 'nV', G.nV);
    fg  = rheome.operators.face_gradient(S.Vertices, S.Faces);
    Pc  = arr.Pos - mean(arr.Pos,1);
    e1 = [1 0 0].'; e2 = [0 1 0].';
    x1 = Pc*e1; x2 = Pc*e2;
    phi = atan2(x2, x1);  r = hypot(x1, x2);
    lam = 0.35*max(r);  k = 2*pi/lam;  w = 2*pi*f0;
    kr  = 2*pi/(0.45*max(r));

    gens = {'travelling','target','rotating','spiral'};
    ph   = { k*x1, kr*r, phi, phi + kr*r };
    nG   = numel(gens);
    % ⚠ THE INITIALISER'S FIELD LIST GOVERNS. Assigning a struct with a field this list
    % lacks errors with "dissimilar structures"; assigning one with FEWER drops them
    % silently. Keep the two in step.
    st = struct('name', {}, 'winding', {}, 'divBulk', {}, 'curlBulk', {}, ...
                'div', {}, 'crl', {}, 'kv', {}, 'core', {});

    fprintf('\n=== rheome.demos.sensor_topology : %s, %d sensors ===\n', arr.Name, arr.nCh);
    for g = 1:nG
        U = cos(ph{g} - w*tv);                       % [nV x nT], the planted generator
        z = i_analytic(U);
        pg = rheome.flow.phasegradient(z(:, jm-2:jm+2), S, 'Rate', fs);
        kf = pg.k(:,:,3);                            % [nF x 3] wavevector per face
        kv = fg.W * kf;                              % to vertices, area-weighted
        kv3 = reshape(kv.', [], 1);                  % [3nV x 1] ambient
        dv  = rheome.differential.divergence(kv3, S, fg);
        cl  = rheome.differential.curl(kv3, S, fg);
        ps  = rheome.detect.phasesingularity(z(:, jm), S);

        % ⚠ MEASURE IN THE BULK, NOT AT THE CORE. |k| diverges as 1/r at a singularity, so
        % div and curl are unbounded there for EVERY rotational pattern -- a peak taken at
        % the core reports the singularity, not the generator. The annulus below also drops
        % the outer edge, where the mesh boundary has no opposite face.
        rc = r;  core = [];
        if ~isempty(ps.charge)
            core = ps.pos(1,:);
            rc = vecnorm(S.Vertices - core, 2, 2);
        end
        R    = max(r);
        bulk = rc > 0.25*R & r < 0.75*R;
        Lk   = median(abs(vecnorm(kv(bulk,:), 2, 2)));       % the local |k| to normalise by
        st(g) = struct('name', gens{g}, ...
            'winding',  sum(abs(ps.charge)), ...
            'divBulk',  median(abs(dv(bulk))) * R / max(Lk, eps), ...
            'curlBulk', median(abs(cl(bulk))) * R / max(Lk, eps), ...
            'div', dv, 'crl', cl, 'kv', kv, 'core', core); %#ok<AGROW>
        fprintf('  %-11s winding %d   div(bulk) %6.3f   |curl|(bulk) %.2e\n', ...
            gens{g}, st(g).winding, st(g).divBulk, st(g).curlBulk);
    end

    % ---------- Fig 1: the wavevector field ----------
    f1 = fbd_fig(doExport, [50 50 1400 380]);
    for g = 1:nG
        subplot(1,nG,g); i_quiv(S, st(g).kv);
        title(gens{g}, 'FontWeight','normal');
    end
    sgtitle('sensor\_topology -- Fig 1: k = \nabla\phi, the vector field a scalar signal does have');
    fbd_save(f1, outDir, 'sensor_topology_1_field');

    % ---------- Fig 2: divergence ----------
    f2 = fbd_fig(doExport, [60 60 1400 380]);
    for g = 1:nG
        subplot(1,nG,g); i_map(x1, x2, st(g).div, true);
        title(sprintf('%s: div k', gens{g}), 'FontWeight','normal');
    end
    sgtitle(['sensor\_topology -- Fig 2: divergence marks where wavefronts are BORN. ' ...
             'Zero for a plane wave; a bright focus for a target and a spiral']);
    fbd_save(f2, outDir, 'sensor_topology_2_divergence');

    % ---------- Fig 3: curl, and the integer that replaces it ----------
    f3 = fbd_fig(doExport, [70 70 1400 380]);
    for g = 1:nG
        subplot(1,nG,g); i_map(x1, x2, st(g).crl, true); hold on;
        if ~isempty(st(g).core)
            c = (st(g).core - mean(S.Vertices,1));
            plot(c(1)*1e3, c(2)*1e3, 'ko', 'MarkerSize', 12, 'LineWidth', 2);
        end
        title(sprintf('%s: curl k  (winding %d)', gens{g}, st(g).winding), 'FontWeight','normal');
    end
    sgtitle(['sensor\_topology -- Fig 3: curl of a gradient is zero EVERYWHERE except at a ' ...
             'singularity. That is the detector, not a defect']);
    fbd_save(f3, outDir, 'sensor_topology_3_curl');

    % ---------- Fig 4: the classification ----------
    f4 = fbd_fig(doExport, [80 80 760 620]);
    hold on; grid on;
    for g = 1:nG
        plot(st(g).winding, st(g).divBulk, 'o', 'MarkerSize', 11, 'LineWidth', 2);
        text(st(g).winding + 0.04, st(g).divBulk, gens{g}, 'FontSize', 11);
    end
    xlim([-0.35 1.6]); xlabel('|winding|  (integer, from the phase circulation)');
    ylabel('bulk |div k| \cdot R / |k|   (dimensionless)');
    title({'two invariants, four families'; ...
           'a spiral is the one carrying BOTH rotation and a source'});
    fbd_save(f4, outDir, 'sensor_topology_4_classify');

    % ---------- checks ----------
    fprintf('\n  self-checks\n');
    gi = @(n) find(strcmp({st.name}, n));
    chk(end+1) = fbd_check('travelling: no winding', st(gi('travelling')).winding, 0, 0.5, '');
    chk(end+1) = fbd_check('target: no winding',     st(gi('target')).winding,     0, 0.5, '');
    chk(end+1) = fbd_check('rotating: one core',     st(gi('rotating')).winding,   1, 0.5, '');
    chk(end+1) = fbd_check('spiral: one core',       st(gi('spiral')).winding,     1, 0.5, '');
    chk(end+1) = fbd_check('travelling: div ~ 0', st(gi('travelling')).divBulk, 0, 0.30, '');
    chk(end+1) = fbd_check('rotating: div ~ 0',   st(gi('rotating')).divBulk,   0, 0.30, '');
    chk(end+1) = fbd_check('target: div > 0',   double(st(gi('target')).divBulk > 0.8), 1, 0.5, '');
    chk(end+1) = fbd_check('spiral: div > 0',   double(st(gi('spiral')).divBulk > 0.8), 1, 0.5, '');
    % The mathematical claim, and the sharpest number in the demo.
    chk(end+1) = fbd_check('curl of a gradient vanishes in the bulk', ...
        max([st.curlBulk]), 0, 1e-6, '');
    chk(end+1) = fbd_check('divergence separates the two pairs', ...
        double(min(st(gi('target')).divBulk, st(gi('spiral')).divBulk) > ...
           3*max(st(gi('travelling')).divBulk, st(gi('rotating')).divBulk)), 1, 0.5, '');
    [ok, ~] = fbd_report(chk);

    out = struct('ok', ok, 'checks', chk, 'stats', st);
end

% ===================== helpers =====================

function z = i_analytic(U)
    U = U - mean(U, 2);  nT = size(U,2);
    F = fft(U, [], 2);  h = zeros(1, nT);  h(1) = 1;
    if mod(nT,2)==0, h(nT/2+1) = 1; h(2:nT/2) = 2; else, h(2:(nT+1)/2) = 2; end
    z = ifft(F .* h, [], 2);
end

function i_quiv(S, kv)
    V = S.Vertices - mean(S.Vertices,1);
    s = 3:4:size(V,1);
    quiver(V(s,1)*1e3, V(s,2)*1e3, kv(s,1), kv(s,2), 1.3, 'LineWidth', 0.9);
    axis equal tight; xlabel('mm');
end

function i_map(x1, x2, v, symm)
    scatter(x1*1e3, x2*1e3, 26, v, 'filled'); axis equal tight; xlabel('mm');
    if symm
        m = prctile(abs(v), 98);  if m <= 0, m = 1; end
        caxis([-m m]);
        n = 128; t = linspace(0,1,n).';
        colormap(gca, [[0.13+0.84*t, 0.29+0.68*t, 0.55+0.42*t]; ...
                       [0.97-0.27*t, 0.97-0.82*t, 0.97-0.81*t]]);
    end
    colorbar;
end

% Author: Diellor Basha, 2026
