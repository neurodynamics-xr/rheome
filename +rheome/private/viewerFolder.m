function vdir = viewerFolder()
% VIEWERFOLDER  The Cortical Flow viewer's folder (dist-viewer) for rheome.view and rheome.snapshot: RHEOME_VIEWER, else
% prefdir/rheome-viewer, installed there from <rheome>/viewer when that ships a newer copy.
    vdir = getenv('RHEOME_VIEWER');
    if ~isempty(vdir), return; end
    vdir = fullfile(prefdir, 'rheome-viewer');
    shipped = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'viewer');
    if isfile(fullfile(shipped, 'viewer.html')) && ~strcmp(stamp(shipped), stamp(vdir))
        if isfolder(vdir), rmdir(vdir, 's'); end
        copyfile(shipped, vdir);
    end
end

function v = stamp(d)
    v = '';
    f = fullfile(d, 'VERSION');
    if isfile(f), v = strtrim(fileread(f)); end
end
