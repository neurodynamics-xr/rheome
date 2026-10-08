function r = outroot()
% LOAD.OUTROOT  Absolute path to the ANALYSIS OUTPUT workspace: results, figures, reports.
%
%   r = rheome.load.outroot()          % $RHEOME_OUT if set, else <repo>/results in a git clone,
%                               % else <userpath>/rheome/results for the installed toolbox
%
% ⭐ TWO ROOTS, TWO KINDS OF FILE. rheome.load.root() is +data: CACHED INPUTS that are regenerable from the
% source recordings (bases, study, envelope modes) and never hand-written. rheome.load.outroot() is where the
% analyses write what they FOUND -- tables, .mat results, figures, rendered reports -- laid out by
% rheome.load.outpath as
%       <root>/<subject>/<analysis>/<file>      a subject's analysis (e.g. a reference subject/alpha_occupancy/)
%       <root>/planted/<analysis>/<file>        subject-free validation on planted or injected data
%       <root>/reports/<file>                   rendered reports; their figures in reports/figures/
% Both are gitignored: the scripts are the record, the outputs are regenerable from them.
%
% ⚠ NOT tempdir, and NOT the current folder. Twenty-odd scripts once wrote to tempdir, which the OS
% wipes on reboot, and every figure and report landed in whatever folder MATLAB happened to be in.
%
% Set RHEOME_OUT (setenv, or the shell) to put the outputs anywhere else without touching a script.
%
% See also: rheome.load.outpath, rheome.load.root
%
% Author: Diellor Basha, 2026

    r = ld_home('RHEOME_OUT', 'results', 'results');
end

% Author: Diellor Basha, 2026
