function out = article(spec)
% REPORT.ARTICLE  Render a filled rheome.report.template as a mini-article: PDF, HTML or Word.
%
%   out = rheome.report.article(spec)
%   spec.Format = "html-file";  out = rheome.report.article(spec)
%
% Builds the standard shape -- title page, contents, abstract, introduction, methods, results
% with numbered figures and legends, conclusion, provenance -- from the spec value, using
% MATLAB Report Generator. Returns the path written.
%
% ⭐ A FIGURE WITHOUT A LEGEND IS AN ERROR, not a warning. An exhibit whose caption says what
% the reader is looking at is the difference between a report and a slide, and the constraint
% is worth enforcing where it is cheap. Set .Legend on every result that sets .Figure.
%
% ⚠ A DIRTY WORKING TREE PUTS A BANNER IN THE DOCUMENT. When spec.Provenance reports
% uncommitted changes the title page carries a warning and the provenance chapter lists the
% files, because a result built from uncommitted code cannot be reproduced from its commit and
% the reader cannot discover that anywhere else.
%
% ⚠ FORMAT CHANGES THE FIGURE SIZING, so check the one you will actually circulate. "pdf" lays
% out to a 6-inch text column; "html-file" is one portable file with images base64-embedded and
% no width constraint; "docx" is for people who will edit it. The same spec gives all three
% (measured on the resolution article: pdf 2.1 s / 537 kB, docx 1.1 s / 505 kB, html 0.8 s /
% 735 kB).
%
% ⚠ IN WORD THE FIGURE AND TABLE NUMBERS ARE FIELDS, NOT TEXT. The caption carries
% "SEQ figure \* ARABIC", and Word computes the number when it updates fields -- which it does
% on open or print, but a reader who never triggers an update can see a stale or zero number.
% Ctrl-A then F9 forces it. The PDF has the numbers baked in and does not have this problem,
% so send the PDF when the numbering has to be right without the reader doing anything.
%
% ⚠ FIGURE HANDLES ARE CAPTURED AT RENDER TIME, so a handle that has been closed or redrawn
% since the result was computed reports whatever it holds NOW. Passing a PNG path that the
% analysis script exported is the reproducible route and is what the example uses.
%
% ⭐ EVERY REPORT IS A SELF-CONTAINED FOLDER. The document goes to spec.OutputPath -- by convention
% rheome.load.outpath("reports/<name>", "<name>"), i.e. results/reports/<name>/<name>.pdf -- and every
% figure it shows is COLLECTED into <that folder>/figures/: a PNG referenced from anywhere else (an
% analysis figure under results/<subject>/...) is copied in, and a figure handle is exported there
% rather than to a throwaway temp file. So a report folder holds the documents and the exact
% images they embed, and can be moved or circulated as one unit.
%
% INPUT:  spec  from rheome.report.template, with the fields filled in
% OUTPUT: out   full path of the document written
%
% See also: rheome.report.template, rheome.report.provenance, mlreportgen.report.Report
%
% Author: Diellor Basha, 2026

    arguments
        spec (1,1) struct
    end
    need = {'Title','Abstract','Introduction','Methods','Results','Conclusion'};
    miss = need(~isfield(spec, need));
    if ~isempty(miss)
        error('report:article:spec', 'spec is missing: %s. Start from rheome.report.template.', ...
              strjoin(miss, ', '));
    end
    fmt = "pdf";  if isfield(spec,'Format') && strlength(spec.Format), fmt = string(spec.Format); end
    outPath = rheome.load.outpath("reports/report", "report");
    if isfield(spec,'OutputPath') && strlength(string(spec.OutputPath))
        outPath = char(spec.OutputPath);
    end
    for k = 1:numel(spec.Results)
        r = spec.Results(k);
        if ~isempty(r.Figure) && (~isfield(r,'Legend') || strlength(string(r.Legend)) == 0)
            error('report:article:legend', ...
                'Result %d has a figure and no legend. Every figure needs one.', k);
        end
    end
    % collect every figure into <report folder>/figures, so the folder is self-contained
    figDir = fullfile(fileparts(outPath), 'figures');
    if ~isfolder(figDir), mkdir(figDir); end
    for k = 1:numel(spec.Results)
        f = spec.Results(k).Figure;
        if isempty(f), continue; end
        if isgraphics(f)
            dst = fullfile(figDir, sprintf('figure_%d.png', k));
            exportgraphics(f, dst, 'Resolution', 150);
        else
            src = i_abs(f);  [~, nm, ex] = fileparts(src);  dst = fullfile(figDir, [nm ex]);
            if ~strcmp(src, dst), copyfile(src, dst); end
        end
        spec.Results(k).Figure = dst;
    end

    import mlreportgen.report.*
    import mlreportgen.dom.*

    % ⚠⚠ PDF CANNOT BE BUILT ON AN exFAT VOLUME, so the document is built under tempdir and
    % the finished file moved in. The chain was measured, not guessed (2026-09-25):
    %   1. An external exFAT volume has has no extended attributes, and macOS sets
    %      com.apple.provenance on every file it writes -- so each one gets an AppleDouble
    %      "._NAME" sibling.
    %   2. Only the PDF path writes an intermediate XSL-FO package DIRECTORY (<name>_FO)
    %      during close(), so that directory fills with "._" siblings, one of which is
    %      "._[Content_Types].xml".
    %   3. The packager enumerates every file in it as a package part. Proven by racing a
    %      deleter against close(): removing the "._" files mid-build changes the error to
    %      "Cannot read. File does not exist", so it was reading them.
    %   4. A part name becomes a URI, and java.net.URI('/._[Content_Types].xml') throws
    %      "Illegal character in path at index 3" -- brackets are URI gen-delims. The real
    %      "[Content_Types].xml" is a reserved OPC name handled separately and is fine.
    %   -> close() fails with "Could not successfully parse URI string", naming neither the
    %      file nor the cause.
    % ⭐ "docx" AND "html-file" BUILD FINE ON exFAT -- measured. docx is a zipped OPC package
    % too, so the earlier guess that "AppleDouble breaks the OPC packager" was wrong; it is
    % specifically the FO intermediate directory that gets scanned. Building everything under
    % tempdir anyway costs nothing and keeps one code path.
    % A failed build also leaves a <name>_FO directory the next open() cannot remove, because
    % MATLAB's recursive delete does not see the "._" entries either -- so that second symptom
    % masks the first. i_clearstale cleans up anything left by an older version of this code.
    i_clearstale(outPath);
    [outDir, outName] = fileparts(outPath);
    if isempty(outDir), outDir = pwd; end
    work = [tempname '_rpt'];
    mkdir(work);
    wc = onCleanup(@() i_rmdir(work));
    R = Report(fullfile(work, outName), char(fmt));
    R.Layout.Landscape = false;
    open(R);
    c = onCleanup(@() i_safeclose(R));

    dirty = isfield(spec,'Provenance') && ~isempty(spec.Provenance) && ...
            spec.Provenance.git.available && ~spec.Provenance.git.clean;

    tp = TitlePage('Title', char(spec.Title));
    if isfield(spec,'Subtitle'), tp.Subtitle = char(spec.Subtitle); end
    if isfield(spec,'Authors'),  tp.Author   = char(strjoin(string(spec.Authors), ', ')); end
    tp.PubDate = char(datetime('today','Format','d MMMM yyyy'));
    add(R, tp);
    if dirty
        p = Paragraph(sprintf(['BUILT FROM AN UNCOMMITTED WORKING TREE (%d modified file(s)). ' ...
            'The commit named in the provenance chapter does NOT reproduce this document.'], ...
            numel(spec.Provenance.git.dirtyFiles)));
        p.Color = 'red';  p.Bold = true;  p.Style = [p.Style {OuterMargin('0in','0in','6pt','12pt')}];
        add(R, p);
    end
    add(R, TableOfContents);

    add(R, i_prose('Abstract', spec.Abstract));
    add(R, i_prose('Introduction', spec.Introduction));
    add(R, i_prose('Methods', spec.Methods));

    ch = Chapter('Title', 'Results');
    for k = 1:numel(spec.Results)
        r = spec.Results(k);
        i_addtext(ch, r.Text);
        if ~isempty(r.Figure)
            fi = FormalImage();
            if isgraphics(r.Figure)
                tmp = [tempname '.png'];
                exportgraphics(r.Figure, tmp, 'Resolution', 150);
                fi.Image = tmp;
            else
                % ⚠ THE DOCUMENT PACKAGER NEEDS AN ABSOLUTE PATH. A relative one fails at
                % close() with "Could not successfully parse URI string", which names neither
                % the image nor the reason.
                fi.Image = char(i_abs(r.Figure));
            end
            fi.Caption = char(r.Legend);
            if fmt == "pdf", fi.ScaleToFit = true; end
            add(ch, fi);
        end
        if isfield(r,'Table') && ~isempty(r.Table)
            bt = BaseTable(MATLABTable(r.Table));
            if isfield(r,'Caption') && strlength(string(r.Caption))
                bt.Title = char(r.Caption);
            end
            add(ch, bt);
        end
    end
    add(R, ch);

    add(R, i_prose('Conclusion', spec.Conclusion));

    if isfield(spec,'Provenance') && ~isempty(spec.Provenance)
        add(R, i_provenance(spec.Provenance));
    end

    close(R);  clear c
    built = R.OutputPath;
    [~, ~, ext] = fileparts(built);
    out = fullfile(outDir, [outName ext]);
    if exist(out, 'file'), delete(out); end
    movefile(built, out, 'f');
end

% ---------------------------------------------------------------------------------------
function ch = i_prose(title, body)
    import mlreportgen.report.Chapter
    ch = Chapter('Title', title);
    i_addtext(ch, body);
end

function i_addtext(container, body)
% Accepts a string, a char, a cellstr or a string array; each element becomes a paragraph.
% ⚠ a placeholder left from rheome.report.template is rendered grey and italic so an unwritten
% section is visible in the document instead of reading as deliberate prose.
    import mlreportgen.dom.*
    if isempty(body), return; end
    s = string(body);
    for i = 1:numel(s)
        if strlength(strtrim(s(i))) == 0, continue; end
        p = Paragraph(char(s(i)));
        if startsWith(strtrim(s(i)), "<") && endsWith(strtrim(s(i)), ">")
            p.Color = '#909090';  p.Italic = true;
        end
        add(container, p);
    end
end

function ch = i_provenance(P)
    import mlreportgen.report.*
    import mlreportgen.dom.*
    ch = Chapter('Title', 'Provenance');
    add(ch, Paragraph(['Collected at run time by rheome.report.provenance; nothing in this chapter ' ...
                       'was typed by hand.']));
    add(ch, Section('Title', 'Run'));
    add(ch, BaseTable(MATLABTable(P.table)));
    if P.git.available && ~P.git.clean
        p = Paragraph(sprintf('Uncommitted at build time (%d):', numel(P.git.dirtyFiles)));
        p.Bold = true;  add(ch, p);
        add(ch, UnorderedList(cellstr(P.git.dirtyFiles)));
    end
    if ~isempty(P.code) && height(P.code) > 0
        add(ch, Section('Title', 'Code'));
        add(ch, BaseTable(MATLABTable(P.code)));
    end
    if ~isempty(P.data) && height(P.data) > 0
        add(ch, Section('Title', 'Data'));
        add(ch, BaseTable(MATLABTable(P.data)));
    end
    add(ch, Section('Title', 'Environment'));
    add(ch, BaseTable(MATLABTable(P.env.toolboxes)));
end

function i_rmdir(p)
    if exist(p, 'dir')
        [ok, ~] = rmdir(p, 's');
        if ~ok && isunix, system(sprintf('rm -rf "%s"', p)); end
    end
end

function i_clearstale(outPath)
    [d, n] = fileparts(outPath);
    if isempty(d), d = pwd; end
    L = dir(fullfile(d, [n '_*']));
    for i = 1:numel(L)
        if ~L(i).isdir, continue; end
        p = fullfile(L(i).folder, L(i).name);
        [ok, ~] = rmdir(p, 's');
        if ~ok && isunix, system(sprintf('rm -rf "%s"', p)); end
    end
end

function a = i_abs(f)
    f = char(f);
    if ~exist(f, 'file')
        error('report:article:figure', 'Figure not found: %s', f);
    end
    d = dir(f);
    a = fullfile(d.folder, d.name);
end

function i_safeclose(R)
    try, close(R); catch, end
end
% Author: Diellor Basha, 2026
