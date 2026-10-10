% Author: Diellor Basha, 2026
classdef tPublicHelp < matlab.unittest.TestCase
% TPUBLICHELP  Rheome is public: its help text and docs name no participant, private script or internal note.
% Scans every .m file under +rheome and doc, every page under doc (generated reference included), README.md and
% Contents.m for participant IDs (sub-0002, sub-MTL0002, sub-PD1201), the private analysis scripts (*_omega, the sphere and
% catalogue validations), the private repository, internal paths (desk/evidence, docs/<date>-*.md), ledger
% request IDs (8 hex digits) and manuscript versions. Examples must stay generic ('mysubject').

properties (Constant)
    Patterns = {
        '[A-Za-z]{2}\w*_omega(?!\w)'                                        % private analysis scripts
        'omega_rerun|sphere_validate|sphere_tilepath_validate|flow_curl_methods|pattern_catalogue_validate|export_data\.m'
        'sub-0\d{3}(?!\d)|sub-[A-Z]{2,}\d+'                                % participant IDs (OMEGA, PREVENT-AD)
        'nxr-cortical-flow'                                                 % the private repository
        'desk/evidence|docs/\d{4}-\d\d-\d\d-|docs/atlas/'                   % internal notes
        '(?<![\w-])(?=[0-9a-f]{0,7}[a-f])(?=[0-9a-f]{0,7}\d)[0-9a-f]{8}(?![\w-])'   % ledger request / task IDs
        '(?<!\w)v18(?!\w)'                                                  % manuscript version
        }
end

methods (Test)
    function publicTextNamesNothingPrivate(tc)
        root = fileparts(fileparts(mfilename('fullpath')));
        files = [dir(fullfile(root, '+rheome', '**', '*.m')); dir(fullfile(root, 'doc', '**', '*.m'))
                 dir(fullfile(root, 'doc', '**', '*.html')); dir(fullfile(root, 'doc', '*.xml'))
                 dir(fullfile(root, 'README.md')); dir(fullfile(root, 'Contents.m'))];
        tc.assertNotEmpty(files);
        bad = strings(0, 1);
        for k = 1:numel(files)
            f = fullfile(files(k).folder, files(k).name);
            txt = regexprep(fileread(f), 'data:[\w/+.-]+;base64,[A-Za-z0-9+/=]+', '');   % embedded figures
            lines = splitlines(string(txt));
            for p = 1:numel(tc.Patterns)
                hit = find(~cellfun(@isempty, regexp(lines, tc.Patterns{p}, 'once')));
                bad = [bad; extractAfter(string(f), strlength(root) + 1) + ":" + hit + ": " + strtrim(regexprep(lines(hit), "^(.{0,120}).*$", "$1"))]; %#ok<AGROW>
            end
        end
        tc.verifyEmpty(bad, "Private details in public help/docs:" + newline + strjoin(unique(bad), newline));
    end

    function patternsCatchTheKnownCases(tc)
        cases = ["see plant_scale_omega.m" "on sub-0002" "PREVENT-AD sub-MTL0002" "sub-PD1201" "nxr-cortical-flow-matlab" ...
                 "desk/evidence/x/plan.md" "See docs/2026-09-22-ingest-notes.md" "(request ea515f28)" "v18 Fig. 8"];
        clean = ["omega_l = l*Omega" "i_omega(fs,nT)" "sub-01_task-noise" "sub-band" "rng(20260905)" "v1.0.0"];
        hits = @(s) any(cellfun(@(p) ~isempty(regexp(s, p, 'once')), tc.Patterns));
        for s = cases, tc.verifyTrue(hits(s), s); end
        for s = clean, tc.verifyFalse(hits(s), s); end
    end
end
end
