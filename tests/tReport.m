classdef tReport < matlab.unittest.TestCase
% The mini-article: the template's shape, the two opinionated rules, and provenance.
%
% ⚠ EVERY TEST BUILDS "html-file", NOT PDF. Same code path through the document packager, and
% 0.3 s against 3 s -- and the output is one greppable file, so the banner and the placeholder
% styling can actually be asserted instead of assumed.
%
% Author: Diellor Basha, 2026

    properties
        tmp; png
    end

    methods (TestClassSetup)
        function build(t)
            % ⚠ NOT under the working directory. Report Generator's packager fails on this
            % volume because macOS writes AppleDouble "._" files into the unpacked package
            % (see rheome.report.article). tempdir is local.
            t.tmp = [tempname '_rpttest'];  mkdir(t.tmp);
            f = figure('Visible','off');  plot(1:10);
            t.png = fullfile(t.tmp, 'fig.png');
            exportgraphics(f, t.png);  close(f);
        end
    end

    methods (TestClassTeardown)
        function clean(t), if exist(t.tmp,'dir'), rmdir(t.tmp,'s'); end, end
    end

    methods (Test)

        function aReportFolderCollectsItsFigures(t)
            % a figure referenced from ANOTHER folder is copied into <report folder>/figures, and a
            % figure HANDLE is exported there, so the report folder is self-contained
            d = fullfile(t.tmp, 'rep_one');  mkdir(d);
            s = i_spec(t);  s.OutputPath = fullfile(d, 'rep_one');
            f = figure('Visible','off');  plot(1:3);
            s.Results = struct('Text', {"a","b"}, 'Figure', {t.png, f}, 'Legend', {"legend a","legend b"}, ...
                               'Table', {[] []}, 'Caption', {"",""});      % cells distribute; strings do not
            rheome.report.article(s);  close(f);
            t.verifyTrue(isfile(fullfile(d, 'figures', 'fig.png')), 'referenced PNG copied in');
            t.verifyTrue(isfile(fullfile(d, 'figures', 'figure_2.png')), 'handle exported in');
            t.verifyTrue(isfile(t.png), 'the original is left where it was');
        end
        function theTemplateHasEverySection(t)
            s = rheome.report.template("A title");
            for f = {'Title','Subtitle','Authors','Abstract','Introduction','Methods', ...
                     'Results','Conclusion','Provenance','Format','OutputPath'}
                t.verifyTrue(isfield(s, f{1}), f{1});
            end
            t.verifyEqual(s.Title, "A title");
            t.verifyTrue(isfield(s.Results, 'Legend'));
        end

        function anIncompleteSpecIsRefused(t)
            s = rheome.report.template();  s = rmfield(s, 'Methods');
            t.verifyError(@() rheome.report.article(s), 'report:article:spec');
        end

        function aFigureWithoutALegendIsAnError(t)
            % ⭐ the rule that makes it a report rather than a slide
            s = i_spec(t);
            s.Results(1).Figure = t.png;
            s.Results(1).Legend = "";
            t.verifyError(@() rheome.report.article(s), 'report:article:legend');
        end

        function aMissingFigureFileIsNamed(t)
            s = i_spec(t);
            s.Results(1).Figure = fullfile(t.tmp, 'nope.png');
            s.Results(1).Legend = "a legend";
            t.verifyError(@() rheome.report.article(s), 'report:article:figure');
        end

        function itBuildsAndTheFigureAndLegendAreInIt(t)
            s = i_spec(t);
            s.Results(1).Figure = t.png;
            s.Results(1).Legend = "A distinctive legend sentence.";
            out = rheome.report.article(s);
            t.verifyEqual(exist(out,'file'), 2);
            txt = fileread(out);
            t.verifyTrue(contains(txt, 'A distinctive legend sentence.'));
            t.verifyTrue(contains(txt, 'Introduction'));
            t.verifyTrue(contains(txt, 'Provenance') || isempty(s.Provenance));
        end

        function aDirtyTreePutsABannerInTheDocument(t)
            % ⚠ the point of the whole provenance chapter: a reader must be told that the
            % commit named in it does not reproduce what they are holding.
            s = i_spec(t);
            s.Provenance = rheome.report.provenance(Hash="none");
            s.Provenance.git.available = true;
            s.Provenance.git.clean = false;
            s.Provenance.git.dirtyFiles = ["a.m"; "b.m"];
            s.Provenance.table = table("working tree", "DIRTY -- 2 uncommitted file(s)", ...
                                       'VariableNames', {'item','value'});
            out = rheome.report.article(s);
            txt = fileread(out);
            t.verifyTrue(contains(txt, 'UNCOMMITTED WORKING TREE'));
            t.verifyTrue(contains(txt, 'a.m'));
        end

        function aCleanTreeGetsNoBanner(t)
            s = i_spec(t);
            s.Provenance = rheome.report.provenance(Hash="none");
            s.Provenance.git.clean = true;
            s.Provenance.git.dirtyFiles = string.empty;
            out = rheome.report.article(s);
            t.verifyFalse(contains(fileread(out), 'UNCOMMITTED WORKING TREE'));
        end

        function provenanceCollectsTheRepositoryAndTheEnvironment(t)
            P = rheome.report.provenance(Hash="none");
            t.verifyTrue(P.git.available);                  % this project IS a git repo
            t.verifyEqual(strlength(P.git.sha), 40);
            t.verifyTrue(islogical(P.git.clean));
            t.verifyFalse(contains(P.git.subject, char(27)));   % ⚠ git's pager escapes, stripped
            t.verifyGreaterThan(height(P.env.toolboxes), 1);
            t.verifyTrue(ismember('working tree', P.table.item));
        end

        function sha256IsStableAndFingerprintIsCheaper(t)
            a = rheome.report.provenance(Data=string(t.png), Hash="sha256");
            b = rheome.report.provenance(Data=string(t.png), Hash="sha256");
            t.verifyEqual(a.data.digest(1), b.data.digest(1));
            t.verifyEqual(strlength(a.data.digest(1)), 64);
            c = rheome.report.provenance(Data=string(t.png), Hash="fingerprint");
            t.verifyTrue(contains(c.data.digest(1), ':'));
            t.verifyTrue(c.data.exists(1));
        end

        function aMissingDataFileIsRecordedNotThrown(t)
            P = rheome.report.provenance(Data="/no/such/file.mat", Hash="sha256");
            t.verifyFalse(P.data.exists(1));
            t.verifyEqual(P.data.bytes(1), 0);
        end
    end
end

function s = i_spec(t)
    s = rheome.report.template("Test article");
    s.Abstract = "Abstract text.";  s.Introduction = "Introduction text.";
    s.Methods = "Methods text.";    s.Conclusion = "Conclusion text.";
    s.Results = struct('Text', "Result text.", 'Figure', [], 'Legend', "", ...
                       'Table', [], 'Caption', "");
    s.Provenance = [];
    s.Format = "html-file";
    s.OutputPath = fullfile(t.tmp, 'article');
end
% Author: Diellor Basha, 2026
