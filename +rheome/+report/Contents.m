% REPORT  A measured result written up as a mini-article, with its provenance collected.
%
% The shape is fixed -- abstract, introduction, methods, results with numbered figures and
% legends, conclusion, provenance -- and the spec that fills it is a VALUE, so one filled
% template renders as PDF, as a single portable HTML file, or as Word without being edited.
%
% Functions:
%   rheome.report.template   - an empty article with every section present and prompted
%   rheome.report.article    - render a filled template (MATLAB Report Generator)
%   rheome.report.provenance - the commit, the dirty files, the code, the data and the environment,
%                       collected at run time rather than typed
%
% ⭐ THE TWO OPINIONATED RULES, both cheap and both load-bearing. A figure without a legend is
% an error. And a document built from an uncommitted working tree carries a red banner saying
% so on its title page, because the commit it names does not reproduce it and a reader has no
% other way to find out.
%
% Needs MATLAB Report Generator; it renders headless: 3 s for a PDF, 0.3 s for single-file HTML
% (R2023b).
%
% ⚠ PDF WILL NOT BUILD ON AN exFAT VOLUME. Only the PDF path writes an intermediate XSL-FO
% package directory, macOS fills it with AppleDouble "._" siblings because exFAT has no
% extended attributes, the packager reads them as package parts, and "._[Content_Types].xml"
% is not a parseable URI. "docx" and "html-file" are unaffected -- rheome.report.article builds
% everything under tempdir and moves the result in, so this is handled, not merely known.
%
% WHERE THE OUTPUT GOES: spec.OutputPath without an extension; rheome.report.article appends the
% format's own and RETURNS the full path. Renders are gitignored -- the spec and the script
% are the record, not the file.
%
% Example: report_resolution_article.m writes up rheome.inverse.resolution as a four-page article.
%
% Author: Diellor Basha, 2026
