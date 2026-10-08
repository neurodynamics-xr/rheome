function demo = cfdemo_study(folder)
%CFDEMO_STUDY  Write a small synthetic Brainstorm study for the documentation examples.
%   DEMO = CFDEMO_STUDY(FOLDER) writes, under FOLDER, the files a Brainstorm protocol holds for
%   one subject: a two-hemisphere cortex (tess_cortex_demo.mat, with a 'Structures' atlas) and a
%   study folder with a channel file, an unconstrained MEG leadfield
%   (headmodel_surf_os_meg_demo.mat), a recording (data_demo.mat) and its noise covariance
%   (noisecov_full.mat). DEMO = CFDEMO_STUDY() writes to fullfile(tempdir, 'cfdemo').
%
%   Everything is synthetic, so the examples run anywhere, without a Brainstorm installation or
%   any participant's data:
%     cortex     two spheres of radius 50 mm (2562 vertices each), one per hemisphere
%     sensors    150 radial magnetometers on a 120 mm helmet
%     leadfield  a current dipole in an infinite homogeneous medium (Biot-Savart), radial field
%     recording  a 10 Hz vortex planted on the left hemisphere, 60 s at 300 Hz, plus white noise
%     noise run  60 s of the same white noise alone (data_noise.mat), the study's empty room
%
%   DEMO is a struct: .cortexFile, .studyDir, .dataName (the arguments rheome.import.surface and
%   rheome.import.study take), .noiseName (the noise run, for rheome.import.study), .vortexCentre [1 x 3]
%   (metres), .vortexVertex (its global vertex index), .sfreq (Hz).
%
%   This is a stand-in for a real study, not a head model: use your own Brainstorm protocol for
%   anything you want to read physiologically.
%
%   Example:
%     demo = cfdemo_study(fullfile(tempdir, 'cfdemo'));
%     assert(isfile(demo.cortexFile) && isfolder(demo.studyDir))
%
%   See also IMPORT.SURFACE, IMPORT.STUDY, IMPORT.BASES, FLOW.CONTEXT.

% Author: Diellor Basha, 2026
if nargin < 1 || isempty(folder), folder = fullfile(tempdir, 'cfdemo'); end
studyDir = fullfile(folder, 'study');
if ~isfolder(studyDir), mkdir(studyDir); end
rng(7);                                                  % same files every time

% -- cortex: two spheres, one per hemisphere, with the Structures atlas that splits them --
[U, Fu] = rheome.geom.icosphere(4);  nh = size(U, 1);  r = 0.05;
cL = [-0.055 0 0.04];  cR = [0.055 0 0.04];
Vertices = [r*U + cL; r*U + cR];
Faces    = [Fu; Fu + nh];
VertNormals = [U; U];
Atlas = struct('Name', 'Structures', 'Scouts', struct( ...
    'Vertices', {1:nh, nh+(1:nh)}, 'Seed', {1, nh+1}, 'Color', {[1 0 0], [0 0 1]}, ...
    'Label', {'Cortex L', 'Cortex R'}, 'Function', {'Mean', 'Mean'}, 'Region', {'LU', 'RU'}, ...
    'Handles', {[], []}));
Comment = 'cfdemo cortex (synthetic)';
cortexFile = fullfile(folder, 'tess_cortex_demo.mat');
save(cortexFile, 'Vertices', 'Faces', 'VertNormals', 'Atlas', 'Comment');

% -- sensors: radial magnetometers on the upper half of a 120 mm sphere --
[S, ~] = rheome.geom.icosphere(3);
S = S(S(:,3) > 0.05, :);  S = S(1:min(150, end), :);  nC = size(S, 1);
loc = 0.12*S + [0 0 0.04];
Channel = struct('Name', compose('MEG%03d', 1:nC), 'Type', 'MEG', 'Comment', '', ...
    'Loc', num2cell(loc', 1), 'Orient', num2cell(S', 1), 'Weight', 1, 'Group', '');   % 1 x nC, as Brainstorm
Comment = 'cfdemo channels (synthetic)';  %#ok<NASGU>
save(fullfile(studyDir, 'channel_demo.mat'), 'Channel', 'Comment');

% -- leadfield: radial B of a current dipole, Biot-Savart, [nC x 3nV] in T/(A.m) --
nV = size(Vertices, 1);  Gain = zeros(nC, 3*nV);  mu = 1e-7;
for c = 1:nC
    d = loc(c,:) - Vertices;  d3 = vecnorm(d, 2, 2).^3;
    Gain(c, :) = reshape((mu * cross(repmat(S(c,:), nV, 1), d, 2) ./ d3 * -1)', 1, []);   % (q x d).n = q.(d x n)
end
HeadModelType = 'surface';  Comment = 'cfdemo Biot-Savart (synthetic)';  %#ok<NASGU>
save(fullfile(studyDir, 'headmodel_surf_os_meg_demo.mat'), 'Gain', 'HeadModelType', 'Comment');

% -- recording: a 10 Hz vortex on the left hemisphere, plus white sensor noise --
sfreq = 300;  Time = (0:60*sfreq-1)/sfreq;
p = [-0.8 -0.3 0.5];  p = p/norm(p);                    % vortex axis: lateral face of the left sphere
[~, iv] = max(U*p');
x = Vertices(1:nh, :) - cL;                              % positions relative to the centre
ang = acos(min(1, U*p'));                                % angle from the vortex axis
w = exp(-(ang/0.45).^2/2);                               % ~22 mm wide on a 50 mm sphere
Jv = w .* cross(repmat(p, nh, 1), x, 2) / r;             % rotation about p: tangential, curl-only
J  = zeros(3*nV, 1);  J(1:3*nh) = reshape(Jv', [], 1);
q  = 20e-9;                                              % 20 nA.m peak per vertex
noise = 5e-15;                                           % 5 fT white noise
F = (Gain*J) * (q*sin(2*pi*10*Time)) + noise*randn(nC, numel(Time));
ChannelFlag = ones(nC, 1);  Comment = 'cfdemo vortex 10 Hz (synthetic)';  nAvg = 1;  %#ok<NASGU>
save(fullfile(studyDir, 'data_demo.mat'), 'F', 'Time', 'ChannelFlag', 'Comment', 'nAvg');
F = noise*randn(nC, numel(Time));  Comment = 'cfdemo noise run (synthetic)';            % an empty-room run
save(fullfile(studyDir, 'data_noise.mat'), 'F', 'Time', 'ChannelFlag', 'Comment', 'nAvg');

NoiseCov = noise^2 * eye(nC);  nSamples = 6000*ones(nC);  Comment = 'cfdemo noise (synthetic)';  %#ok<NASGU>
save(fullfile(studyDir, 'noisecov_full.mat'), 'NoiseCov', 'nSamples', 'Comment');

demo = struct('cortexFile', cortexFile, 'studyDir', studyDir, 'dataName', 'data_demo.mat', 'noiseName', 'data_noise.mat', ...
    'vortexCentre', cL + r*p, 'vortexVertex', iv, 'sfreq', sfreq);
end
% Author: Diellor Basha, 2026
