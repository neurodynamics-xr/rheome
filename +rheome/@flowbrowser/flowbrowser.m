classdef flowbrowser < handle
% FLOWBROWSER  Interactive browser for the band- and scale-resolved cortical flow.
%
%   B  = rheome.flowpage.prepare('sub01');    % ~40 s, once
%   app = rheome.flowbrowser(B);                      % opens the window
%   app = rheome.flowbrowser(B, 'Band', 'beta', 'PageIndex', 5);
%
% Select a band, advance through pages and scrub through time, and watch the vorticity
% decomposed across spatial scales while the lambda spectrum shows how globally or locally
% the energy is distributed.
%
% ⭐ WHY SCRUBBING IS INSTANT. A page's coefficients are computed once (CWT + one GEMM,
% ~2-3 s); every frame after that is a GEMV per scale, Phi*(h_g .* c(:,t)), ~33 MFLOP. The
% surfaces are built once and only their FaceVertexCData is rewritten, so nothing re-renders.
%
% ⭐ GLOBAL vs LOCAL IS SPECTRAL. Energy in low lambda is globally distributed structure;
% energy in high lambda is spatially confined structure. A vortex occurring while the low
% modes are occupied is a different object from the same vortex while activation is confined,
% and the Spectrum tab is where that distinction is read.
%
% ⚠ THE MARGIN IS SHADED, NOT HIDDEN. Frames outside the page's core belong to the
% neighbouring pages too and carry the wavelet cone; the Timeline shades them and the title
% says so. Scrubbing there is allowed -- believing it uncritically is not.
%
% TABS
%   Cortex     the bank total and one surface per spatial scale, each labelled with its
%              wavelength and its share of the frame's energy, the peak one marked.
%              Set ShowDominant = true for an extra panel colouring each VERTEX by which
%              scale wins there (hue) and how much energy is there (brightness) -- off by
%              default because nine panels crowd the tab; the code stays live.
%              Signed (CCW +/CW -) or Magnitude (the envelope) -- see @flowpage/map.
%   Spectrum   energy per LBO mode and per scale at the current frame, centroid marked.
%   Timeline   the cortical SCALOGRAM (scale x time, column-normalised), band power with the
%              lambda-centroid, and the global index -- all on one time axis.
%
% KEYBOARD  left/right step a frame (shift: ten), up/down change page, space plays,
%           m toggles signed/magnitude.
%
% Adding a view is a new tab function of (flowpage, frame); nothing else changes.
%
% See also: rheome.flowpage, rheome.flowpage.prepare, rheome.show.surface
%
% Author: Diellor Basha, 2026

    properties
        Bundle
        Page                                   % current flowpage
        Frame       (1,1) double = 1
        BandName    (1,:) char = 'alpha'
        PageIndex   (1,1) double = 1
        MapMode     (1,:) char = 'signed'
        NumSensorTraces (1,1) double = 6
        ShowDominant    (1,1) logical = false   % the dominant-scale map: off, it crowds the tab
        Bands
    end

    % Read-only handles: public get so the app can be inspected from a script or a test,
    % private set so nothing outside can swap them out from under the callbacks.
    properties (SetAccess = private)
        Fig                                    % the uifigure
        Tabs                                   % the uitabgroup
        FrameLbl                               % frame / time / margin readout
        SensorAx                               % the sensor-trace strip on the Cortex tab
        SensorCursor = []                      % its time cursor
        ScalogramData = []                     % [nScale x nFrames], cached per page
    end

    properties (Access = private)
        Ctl
        BandDrop; PageSpin; FrameSlider; ModeSwitch; PlayBtn; StatusLbl
        CortexAx = []; CortexPatch = []; CortexTitle = []
        SpecAxMode; SpecAxScale
        TimeAxPower; TimeAxScale; ScalogramAx; TimeCursors = []
        ScaleColors = []; ScaleOrder = []
        Playing (1,1) logical = false
    end

    methods
        function app = flowbrowser(B, varargin)
            p = inputParser;
            p.addParameter('Band',      'alpha');
            p.addParameter('PageIndex', 1);
            p.parse(varargin{:});

            app.Bundle = B;
            app.Bands  = struct('delta',[1 4], 'theta',[4 8], 'alpha',[8 13], ...
                                'beta',[13 30], 'gamma',[30 60]);
            app.BandName  = p.Results.Band;
            app.PageIndex = p.Results.PageIndex;

            app.buildUI();
            app.loadPage();
        end

        function delete(app)
            app.Playing = false;
            if ~isempty(app.Fig) && isvalid(app.Fig), delete(app.Fig); end
        end
    end

    % The callbacks are also the PROGRAMMATIC API -- drive the browser from a script, or
    % from a test, without synthesising UI events.
    methods
        onBand(app, name)
        onPage(app, k)
        onFrame(app, i)
        onMode(app, m)
        onPlay(app, on)
        onKey(app, e)
        setFrame(app, i)
        stepFrame(app, d)
        loadPage(app)
        refresh(app)
    end

    methods (Static)
        c = coolwarm(n)
    end

    methods (Access = private)
        buildUI(app)
        buildCortexAxes(app)
        drawCortex(app)
        drawSpectrum(app)
        drawTimeline(app)
        drawScalogram(app)
        drawDominant(app, Z)
        drawSensors(app)
        updateFrameLabel(app)
        shadeMargin(app, ax)
        updateCursors(app)
        setStatus(app, s)
    end
end

% Author: Diellor Basha, 2026
