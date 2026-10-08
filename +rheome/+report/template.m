function spec = template(title)
% REPORT.TEMPLATE  An empty mini-article, with every section present and prompted.
%
%   spec = rheome.report.template()
%   spec = rheome.report.template("What the array resolves on the cortex")
%   spec.Abstract = "..."; spec.Results(1).Figure = "fig.png"; rheome.report.article(spec)
%
% ⭐ THE TEMPLATE IS A VALUE, not a file to copy. Fill the fields, hand it to rheome.report.article,
% and the same spec can be rebuilt as PDF, single-file HTML or Word without editing anything.
% Every field starts as a PROMPT: whatever is left unedited appears in the document in a
% placeholder style, so an unwritten section is visible rather than silently missing.
%
% FIELDS
%   Title, Subtitle, Authors        strings; Authors may be a string array
%   Abstract                        what was done and what came out, in a paragraph
%   Introduction                    why the question was asked and what it bears on
%   Methods                         how it was measured, in enough detail to repeat
%   Results                         [1 x n] struct array, one entry per exhibit:
%                                     .Text    the prose for this result
%                                     .Figure  path to an image, or a figure handle, or []
%                                     .Legend  ⚠ REQUIRED whenever .Figure is set
%                                     .Table   a MATLAB table, or []
%                                     .Caption caption for .Table
%   Conclusion                      what it means and what it forbids
%   Provenance                      rheome.report.provenance(); [] omits the chapter
%   Format                          "pdf" (default) | "html-file" | "docx"
%   OutputPath                      where to write, without the extension
%
% See also: rheome.report.article, rheome.report.provenance
%
% Author: Diellor Basha, 2026

    arguments
        title (1,1) string = "Untitled result"
    end

    spec = struct();
    spec.Title    = title;
    spec.Subtitle = "<subtitle: the finding in one line>";
    spec.Authors  = "Diellor Basha";
    spec.Abstract = "<abstract: the question, the measurement, the number, the consequence>";
    spec.Introduction = "<introduction: why this was asked, and what rests on the answer>";
    spec.Methods  = "<methods: what was computed, on what data, by which function>";
    spec.Results  = struct('Text', "<result: what the exhibit shows>", ...
                           'Figure', [], 'Legend', "", 'Table', [], 'Caption', "");
    spec.Conclusion = "<conclusion: what this licenses, and what it rules out>";
    spec.Provenance = [];
    spec.Format     = "pdf";
    spec.OutputPath = fullfile(pwd, matlab.lang.makeValidName(title));
end
% Author: Diellor Basha, 2026
