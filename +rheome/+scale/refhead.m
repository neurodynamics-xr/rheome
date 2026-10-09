function file = refhead(name, file)
% SCALE.REFHEAD  A cached head, without its recording, as the reference head of measure_trackfactorial.
%
%   file = rheome.scale.refhead('subject01', '/path/subject01.tar.gz')   % the name must end in .tar.gz
%
% Bundles what rheome.scale.sensors reads -- bases.mat, surface.mat (with Reg.Sphere) and study.mat
% (channels, gain, noise covariance; rec.F emptied) -- from rheome.load.root()/<name>. Unpacked into
% <RHEOME_DATA>/<name> with RHEOME_REFHEAD=<name> set, it is MS1 G8's "ref" head level (nsp cf-track,
% param ref_head; cfscale_subject.sh unpacks it). The tarball is named <name>.tar.gz by convention.
%
% See also: rheome.scale.measure_trackfactorial, rheome.scale.sensors
%
% Author: Diellor Basha, 2026

    d = fullfile(rheome.load.root(), name);
    Sf = rheome.load.surface(name);
    if ~isfield(Sf, 'Sphere') || isempty(Sf.Sphere), error('scale:refhead:noSphere', '%s has no Reg.Sphere.', name); end
    t = tempname;  mkdir(fullfile(t, name));  c = onCleanup(@() rmdir(t, 's'));
    copyfile(fullfile(d, 'bases.mat'), fullfile(t, name));  copyfile(fullfile(d, 'surface.mat'), fullfile(t, name));
    st = load(fullfile(d, 'study.mat'));  st.rec.F = zeros(size(st.rec.F, 1), 0);   % the head, not the recording
    save(fullfile(t, name, 'study.mat'), '-struct', 'st', '-v7.3');
    tar(file, '*.mat', fullfile(t, name));                    % gzipped: the name ends in .gz
end

% Author: Diellor Basha, 2026
