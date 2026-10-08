function p = outpath(analysis, file, subject)
% LOAD.OUTPATH  Where an analysis output goes: <outroot>/[<subject>/]<analysis>/<file>, folder created.
%
%   p = rheome.load.outpath("alpha_occupancy", "alpha_occupancy.mat", "sub01")
%       -> <outroot>/a reference subject/alpha_occupancy/alpha_occupancy.mat
%   p = rheome.load.outpath("plant_tilepath", "plant_tilepath_seed13.csv", "planted")
%   p = rheome.load.outpath("reports", "grouptrack_article")                  % -> <outroot>/reports/grouptrack_article
%   p = rheome.load.outpath("reports/figures", "")                           % the folder itself
%
% The SAME call is used to write a file and, in any other script, to read it, so a producer and its
% consumers cannot disagree about where it lives. The folder is created on every call (cheap, and
% it means a reader never fails on a missing directory -- only on a missing file, which is the error
% worth seeing).
%
% INPUTS
%   analysis  the producing analysis (the script's stem without "_omega"); may contain "/"
%   file      the file name, or "" for the folder
%   subject   the dataset name, "planted" for subject-free validation, or omitted for neither
%
% See also: rheome.load.outroot, rheome.load.root
%
% Author: Diellor Basha, 2026

    arguments
        analysis (1,1) string
        file (1,1) string = ""
        subject (1,1) string = ""
    end
    parts = {rheome.load.outroot()};
    if strlength(subject) > 0, parts{end+1} = char(subject); end
    parts = [parts, strsplit(char(analysis), '/')];
    d = fullfile(parts{:});
    if ~isfolder(d), mkdir(d); end
    if strlength(file) > 0, p = fullfile(d, char(file)); else, p = d; end
end

% Author: Diellor Basha, 2026
