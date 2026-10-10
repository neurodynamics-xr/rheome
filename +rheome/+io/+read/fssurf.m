function S = fssurf(file)
% IO.READ.FSSURF  Read a FreeSurfer binary triangle surface (lh.white, lh.sphere, ...) into the house struct.
%
%   S = rheome.io.read.fssurf(fullfile(subjDir, 'surf', 'lh.sphere'))
%
% The FreeSurfer triangle format: magic 0xFFFFFE (3 bytes), two newline-terminated comment lines, then
% big-endian int32 nV and nF, float32 vertex coordinates [nV x 3] in MILLIMETRES, int32 faces [nF x 3],
% 0-based. ⚠ Coordinates are returned in METRES, like every other surface in the toolbox.
% Vertex i of lh.sphere, lh.white, lh.pial and lh.inflated of one subject is the same point, so a field
% computed on one is drawn on any other without resampling.
%
% OUTPUT (struct S)   .Vertices [nV x 3] m   .Faces [nF x 3] 1-based   .nV .nF   .Comment   .SurfaceFile
%
% See also: rheome.io.read.surface, rheome.operators.face_gradient
%
% Author: Diellor Basha, 2026

    fid = fopen(file, 'r', 'ieee-be');
    if fid < 0, error('io:read:fssurf:notFound', 'Cannot open %s', file); end
    c = onCleanup(@() fclose(fid));
    magic = fread(fid, 3, 'uchar').';
    if ~isequal(magic, [255 255 254])
        error('io:read:fssurf:format', '%s is not a FreeSurfer triangle surface (magic %s)', file, mat2str(magic));
    end
    comment = strtrim(fgetl(fid));  fgetl(fid);           % the creator line, then an empty line
    n = fread(fid, 2, 'int32');
    V = fread(fid, [3 n(1)], 'float32').';
    F = fread(fid, [3 n(2)], 'int32').' + 1;
    if size(F,1) ~= n(2) || any(F(:) < 1 | F(:) > n(1))
        error('io:read:fssurf:format', '%s is truncated or its faces index outside the vertices', file);
    end
    S = struct('Vertices', 1e-3*double(V), 'Faces', double(F), 'nV', n(1), 'nF', n(2), ...
               'Comment', comment, 'SurfaceFile', char(file));
end

% Author: Diellor Basha, 2026
