% VIZ  Rheome's visualization package: every cortical visual is drawn by the Cortical Flow viewer.
%
% Diellor, 2026-10-14: the viewer (WebGPU; shipped in <rheome>/viewer, versioned with Rheome's releases) is Rheome's
% renderer of the cortex. MATLAB graphics (rheome.show) are the fallback -- no WebGPU -- tests and quick looks;
% nothing in a figure pipeline calls them. Everything reads the lab's figure tokens.
%
% THE PRIMARY API
%   rheome.view       interactive view: Layers= (field, form, representation, look), Labellings=, Points=, playback,
%                     .set, .screenshot (PNG + its marks sidecar)
%   rheome.snapshot   print: the figure system's four views, strips (image series over time), the geodesic scale bar,
%                     colourbars by quantity, SVG/PDF/PNG, each render's marks sidecar
%   rheome.io.write.bundle   the data the viewer reads (rheome.bundle@1: surfaces, fields, labellings, point sets)
%
% WHAT THE VIEWER IS GIVEN (this package)
%   rheome.viz.tiles   a tile level of a cortical pyramid (rheome.geom.tree) as a labelling: tiles by level
%   rheome.viz.points  critical points (rheome.flow.critical, rheome.detect.criticalPoints) as point sets, one a type
%   rheome.viz.marks   a render's marks sidecar (<render>.marks.json, nxr.marks@1): every overlay with its key, what
%                      figure_check.py E1 reads (Chekhov's gun)
%
% EVERY CORTICAL VISUAL, AND WHO DRAWS IT (the viewer's README: nxr-cortical-flow editor/src/viewer/README.md)
%   a real field, a vector field's |v| n.v div curl     viewer: a colour / isolines / halfMax layer
%   tangent field (arrows), full field or per tile      viewer: an arrows layer (form x | xt; level)
%   Hodge potentials phi, psi                            bundle 'Potentials' (solved in MATLAB), viewer: colour
%   surfaces (folded, inflated, sphere), hemispheres     viewer: 'surface', 'hemisphere', the folded one primary
%   lighting and silhouette                              viewer (its renderer)
%   four standard views, layouts, theme, key band        rheome.snapshot / the viewer's 'figure' event
%   image series / strips over time                      rheome.snapshot (Time= a vector)
%   geodesic scale bar on the folded surface             rheome.snapshot, viewer 'scaleBar'
%   tiles by level (fill)                                rheome.viz.tiles + viewer 'labelling'
%   critical points                                      rheome.viz.points + viewer 'points'
%   keyed marks (marks sidecar)                          viewer screenshot/figure 'marks' + rheome.viz.marks
%   multiscale insets (a region re-rendered, outlined)   rheome.snapshot Inset= (PR #79)
%   GAPS, filed with Product (the viewer): curvature shading under the data; tile OUTLINES with major/minor weights;
%   a geodesic-path layer (paths, trajectories, arrowhead); selection fills / outlines (a tile, a mask, a support);
%   overlays returned in pixel coordinates (marks sidecar with positions); isolines dashed below zero; the tokens'
%   critical-point glyphs (+ - o). Until each lands, a figure that needs it says so; it is not drawn in rheome.show.
%
% NOT CORTICAL (stay in rheome.show, by the figure tokens): spectra, scalograms, time-tracks, phasor glyphs, the
% level stack, the resolution ruler, ridges -- plots, not surfaces.
%
% FALLBACK, TESTS AND QUICK LOOKS: rheome.show.surface, rheome.show.vectors, rheome.show.filmstrip,
% rheome.show.vectorstrip, rheome.show.gif, rheome.show.frame, rheome.show.critical (MATLAB graphics, no WebGPU).
%
% Author: Diellor Basha, 2026
