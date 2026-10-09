function out = bstopticalflow(X, S, fs, opts)
% FLOW.BSTOPTICALFLOW  Brainstorm's cortical optical flow (bst_opticalflow, Lefevre & Baillet 2008), called as Brainstorm calls it.
%
%   out = rheome.flow.bstopticalflow(X, S, fs)                        % Horn-Schunck 0.01, Brainstorm's default
%   out = rheome.flow.bstopticalflow(X, S, fs, HornSchunck=0.1, BrainstormDir="/path/to/brainstorm3")
%
% The comparator of MS1 G16 (estimator ii). Brainstorm's optical-flow panel (panel_opticalflow.m) reads
% the source map, takes abs(ImageGridAmp) and passes it to bst_opticalflow with the Horn-Schunck weight
% typed in its text field (default 0.01). This does the same and nothing more: F = abs(X), the frames
% 2..nT against their predecessors, and Brainstorm's own code computes the flow. bst_opticalflow
% normalises F by its maximum, so the field is per frame in metres (the mesh's units) and is scaled
% here by fs to m/s.
%
% ⚠ BRAINSTORM IS GPL AND RHEOME IS MIT, so the function is CALLED from a Brainstorm checkout, never
% copied: BrainstormDir (default env RHEOME_BRAINSTORM) is put on the path (toolbox/math and
% toolbox/core) when bst_opticalflow is not already there. Outside a Brainstorm session, GlobalData is
% set to Brainstorm's server mode (Program.GuiLevel = -1), in which bst_progress returns at once; a
% running Brainstorm session is left as it is.
%
% INPUTS   X [nV x nT] source time series (signed; abs is taken), S surface (.Vertices .Faces
%          .VertNormals), fs Hz.
% OUTPUT   out.velocity [nV x 3 x nT-1] m/s   out.poincare [nF x nT-1]   out.errorData out.errorReg
%          out.hornSchunck   out.brainstorm (the bst_opticalflow file used)
%
% See also: rheome.dynamics.opticalflow_scalar (our port, not this), rheome.flow.phaseregression,
%           rheome.scale.measure_catalogue
%
% Author: Diellor Basha, 2026

    arguments
        X {mustBeNumeric}
        S struct
        fs (1,1) double {mustBePositive}
        opts.HornSchunck (1,1) double {mustBePositive} = 0.01
        opts.BrainstormDir string = string(getenv('RHEOME_BRAINSTORM'))
    end
    if exist('bst_opticalflow', 'file') ~= 2
        if strlength(opts.BrainstormDir) == 0 || ~isfolder(fullfile(opts.BrainstormDir, 'toolbox', 'math'))
            error('flow:bstopticalflow:path', ['bst_opticalflow is not on the path: set BrainstormDir (or ' ...
                  'RHEOME_BRAINSTORM) to a Brainstorm checkout']);
        end
        addpath(fullfile(opts.BrainstormDir, 'toolbox', 'math'), fullfile(opts.BrainstormDir, 'toolbox', 'core'));
    end
    global GlobalData %#ok<GVMIS> Brainstorm's own global, read by bst_progress
    if isempty(GlobalData), GlobalData.Program.GuiLevel = -1; end
    FV = struct('Faces', double(S.Faces), 'Vertices', double(S.Vertices), 'VertNormals', double(S.VertNormals));
    nT = size(X, 2);  Time = (0:nT-1) / fs;
    [flow, ~, eD, eR, pc] = bst_opticalflow(abs(double(X)), FV, Time, Time(2), Time(end), opts.HornSchunck);
    out = struct('velocity', flow * fs, 'poincare', pc, 'errorData', eD, 'errorReg', eR, ...
                 'hornSchunck', opts.HornSchunck, 'brainstorm', which('bst_opticalflow'));
end

% Author: Diellor Basha, 2026
