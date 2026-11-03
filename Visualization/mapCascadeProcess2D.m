function mapCascadeProcess2D(Boundary, System, CascadeTrace, Mode, options)
% INTRODUCTION
%   Animate the cascading-failure process of multi-infrastructure systems and
%   (optionally) export the animation as a single GIF.
%   - This function is UI-simple: users only provide Boundary/System/Trace/Mode and a
%     minimal options struct; all other styling options are forwarded to mapMultipleCISLayout2D.
%   - Internally, it creates and reuses ONE figure/axes; each frame is redrawn onto
%     the same axes to produce smooth animation. Axis limits are locked after the
%     first frame to avoid visual zoom in/out between frames.

% INPUT
%   Boundary      : struct for map boundary (same as mapMultipleCISLayout2D).
%                   Common fields: .X, .Y (vectors) for the study boundary.
%   System        : struct with the systems to draw (same schema as mapMultipleCISLayout2D).
%                   Typical fields:
%                     .PowerSystem / .GasSystem / .WaterSystem (each with .Node/.Edge arrays)
%                     .PowerToGas / .GasToPower / .PowerToWater / .WaterToPower (interdependency)
%   CascadeTrace  : 1-by-R struct array, each element describes one frame in time order.
%                   Required fields depend on Mode:
%                     Mode = 'Flow':
%                       (per system)  <Sys>EdgeState, <Sys>DemNode, <Sys>GenNode
%                       (links)        PowerToGas, GasToPower, PowerToWater, WaterToPower (optional)
%                       Flow colors use each system's peak across all frames;
%                       node sizes are normalized per system and frame.
%                     Mode = 'Connectivity':
%                       (per system)  <Sys>NodeState
%                       (links)        PowerToGas, GasToPower, PowerToWater, WaterToPower (optional)
%                       MetricLimits defaults to [0,1], and ValueLabel defaults
%                       to 'Component functionality'.
%                   Here <Sys> is one of {'Power','Gas','Water'}.
%   Mode          : char, REQUIRED. Either 'Flow' or 'Connectivity'.
%                   - 'Flow' calls mapMultipleCISLayout2D with Mode='Flow' and custom flows.
%                   - 'Connectivity' calls mapMultipleCISLayout2D with Mode='Metric' and per-node metrics.
%   options       : struct (minimal options consumed here; all others forwarded to the renderer)
%                   Consumed here:
%                     - FigureColor : (default 'white') figure background color
%                     - FramePause  : (default 1) finite nonnegative seconds per frame;
%                                     also used as GIF DelayTime
%                     - GifFilename : (default 'CascadingProcess2D') base filename without extension;
%                                     if empty '', no file is written. The function saves '<GifFilename>.gif'.
%                   Other fields are passed through to mapMultipleCISLayout2D,
%                   except AxesOpt is replaced by the internal animation axes and
%                   FlowSource is forced to 'custom' in Flow mode.
% 
% OUTPUT
%   (none) This function renders an on-screen animation and optionally writes a GIF file.

if nargin < 4 || isempty(Mode)
    error('Mode is required and must be "Flow" or "Connectivity".');
end
if nargin < 5 || isempty(options), options = struct(); end
Mode = string(Mode);
if ~isscalar(Mode) || ~any(Mode == ["Flow","Connectivity"])
    error('Mode must be "Flow" or "Connectivity".');
end

R = numel(CascadeTrace);
if R==0
    warning('CascadeTrace is empty. Nothing to animate.');
    return;
end

FigureColor = gopt(options,'FigureColor','white');
FramePause = gopt(options,'FramePause',1);
if ~isnumeric(FramePause) || ~isscalar(FramePause) || ~isreal(FramePause) || ...
        ~isfinite(FramePause) || FramePause < 0
    error('options.FramePause must be a finite nonnegative numeric scalar.');
end
if isfield(options,'GifFilename')
    GifFilename = options.GifFilename;
else
    GifFilename = 'CascadingProcess2D';
end
if isstring(GifFilename) && isscalar(GifFilename)
    GifFilename = char(GifFilename);
elseif ~isempty(GifFilename) && (~ischar(GifFilename) || ~isrow(GifFilename))
    error('options.GifFilename must be a character vector, string scalar, or empty string.');
end

hfig = figure('Color', FigureColor, 'Name','Cascade Animation');
hax  = axes('Parent', hfig);

delay = max(0.01, FramePause);
sysList = {'Power','Gas','Water'};
ipsList = {'PowerToGas','GasToPower','PowerToWater','WaterToPower'};
hasSystem = cellfun(@(name) isfield(System,[name 'System']) && ...
    ~isempty(System.([name 'System'])), sysList);
hasIps = cellfun(@(name) isfield(System,name) && ~isempty(System.(name)), ipsList);

if strcmp(Mode,'Flow')
    PaintMode = 'Flow';
    options.FlowSource = 'custom';
    for s = 1:numel(sysList)
        if ~hasSystem(s), continue; end
        field = [sysList{s} 'EdgeState'];
        flowMax = 0;
        for r = 1:R
            values = abs(CascadeTrace(r).(field));
            values = values(isfinite(values));
            if ~isempty(values), flowMax = max(flowMax, max(values)); end
        end
        if flowMax <= 0, flowMax = 1; end
        options.([sysList{s} 'FlowMax']) = flowMax;
    end
else
    PaintMode = 'Metric';
    if ~isfield(options,'MetricLimits') || isempty(options.MetricLimits)
        options.MetricLimits = [0,1];
    end
    if ~isfield(options,'ValueLabel') || isempty(options.ValueLabel)
        options.ValueLabel = 'Component functionality';
    end
end
options.AxesOpt = hax;

for r = 1:R
    opts = options;
    if strcmp(Mode,'Flow')
        for s = 1:numel(sysList)
            if hasSystem(s)
                opts.([sysList{s} 'EdgeFlow']) = CascadeTrace(r).([sysList{s} 'EdgeState']);
                opts.([sysList{s} 'NodeDem']) = CascadeTrace(r).([sysList{s} 'DemNode']);
                opts.([sysList{s} 'NodeGen']) = CascadeTrace(r).([sysList{s} 'GenNode']);
            end
        end
        for i = 1:numel(ipsList)
            if hasIps(i) && isfield(CascadeTrace(r),ipsList{i})
                opts.([ipsList{i} 'Flow']) = interdependencyFlow( ...
                    System.(ipsList{i}), ipsList{i}, CascadeTrace(r).(ipsList{i}));
            end
        end
    else
        for s = 1:numel(sysList)
            if hasSystem(s)
                opts.([sysList{s} 'NodeMetric']) = CascadeTrace(r).([sysList{s} 'NodeState']);
            end
        end
        for i = 1:numel(ipsList)
            if hasIps(i) && isfield(CascadeTrace(r),ipsList{i})
                opts.([ipsList{i} 'Metric']) = double(CascadeTrace(r).(ipsList{i}));
            end
        end
    end

    cla(hax);
    mapMultipleCISLayout2D(Boundary, System, PaintMode, opts);

    if r == 1
        axis(hax,'tight');
        xlim0 = xlim(hax); ylim0 = ylim(hax);
        set(hax,'XLim',xlim0,'YLim',ylim0);
    else
        set(hax,'XLim',xlim0,'YLim',ylim0);
    end

    drawnow;

    if ~isempty(GifFilename)
        [A,map] = rgb2ind(frame2im(getframe(hfig)),256);
        if r==1
            imwrite(A,map,[GifFilename '.gif'],'LoopCount',inf,'DelayTime',delay);
        else
            imwrite(A,map,[GifFilename '.gif'],'WriteMode','append','DelayTime',delay);
        end
    end
    pause(FramePause);
end
end

function val = gopt(s, name, default)
if isfield(s,name) && ~isempty(s.(name)), val = s.(name); else, val = default; end
end

function flow = interdependencyFlow(linkTable, field, state)
if strcmp(field,'GasToPower'), demandCol = 5; else, demandCol = 3; end
flow = state(:).*linkTable(:,demandCol);
end
