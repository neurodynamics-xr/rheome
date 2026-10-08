function G = gfbCaptureGolden(outFile)
% GFBCAPTUREGOLDEN  Snapshot rheome.filters.frame's output for the shim regression.
%
%   G = gfbCaptureGolden()                 % just return it
%   G = gfbCaptureGolden(outFile)          % and save to a .mat
%
% ⚠ RUN THIS ONLY AGAINST THE PRE-SHIM rheome.filters.frame. tests/gfbGoldenFrame.mat was
% captured at commit 4f18db7..., BEFORE rheome.filters.frame was converted to delegate to
% graphfilterbank, and is the reference proving the conversion changed nothing.
% Re-running it after the shim lands would capture the shim's own output and destroy
% the very thing it is meant to check.
%
% Numeric fields only -- function handles are deliberately reduced to evaluated gains,
% so the .mat carries no captured workspace.
%
% Author: Diellor Basha, 2026

    fx  = gfbFixture();
    lam = fx.Lambda;

    cases = {};
    cases{end+1} = {'mexhat',   7,  lam,        {}};
    cases{end+1} = {'mexhat',   [], lam,        {}};
    cases{end+1} = {'mexhat',   [], lam,        {'sigmaMax', 0.06}};
    cases{end+1} = {'mexhat',   5,  lam,        {'tmin','usable'}};
    cases{end+1} = {'itersine', 8,  lam,        {}};
    cases{end+1} = {'heat',     6,  lam,        {}};
    cases{end+1} = {'mexhat',   6,  [200 9318], {}};
    cases{end+1} = {'mexhat',   6,  [0   9318], {}};
    cases{end+1} = {'mexhat',   [], lam,        {'perOctave', 5}};
    cases{end+1} = {'mexhat',   [], lam,        {'lpfactor', 50}};

    G = struct('family',{},'Nf',{},'lrange',{},'opts',{},'M',{},'H',{}, ...
               'Lrange',{},'Centers',{},'t',{},'Sigma',{},'Gamma',{},'SigmaTrue',{}, ...
               'MassLost',{},'Usable',{},'SigmaFloor',{},'SigmaMax',{},'LminEff',{}, ...
               'PerOctave',{},'NfField',{});
    for i = 1:numel(cases)
        c = cases{i};
        f = rheome.filters.frame(c{1}, c{2}, c{3}, c{4}{:}, 'warn', false);
        G(i).family = c{1};  G(i).Nf = c{2};  G(i).lrange = c{3};  G(i).opts = c{4};
        G(i).M       = numel(f.g);
        G(i).H       = rheome.filters.frame_gains(f, lam);
        G(i).NfField = f.Nf;
        for fld = {'Lrange','Centers','t','Sigma','Gamma','SigmaTrue','MassLost', ...
                   'Usable','SigmaFloor','SigmaMax','LminEff','PerOctave'}
            G(i).(fld{1}) = f.(fld{1});
        end
    end

    if nargin >= 1 && ~isempty(outFile)
        save(outFile, 'G', '-v7');
        fprintf('gfbCaptureGolden: %d cases -> %s\n', numel(G), outFile);
    end
end

% Author: Diellor Basha, 2026
