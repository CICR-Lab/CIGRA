function mapCascadeProcess3D(Boundary, System, CascadeTrace, Mode, options)
% INTRODUCTION
%   Animate the cascading-failure process of multi-infrastructure systems and
%   (optionally) export the animation as a single GIF.
%   ? This function is UI-simple: users only provide Boundary/System/Trace/Mode and a
%     minimal options struct; all other styling options are forwarded to mapMultipleCISLayout3D.
%   ? Internally, it creates and reuses ONE figure/axes; each frame is redrawn onto
%     the same axes to produce smooth animation. Axis limits are locked after the
%     first frame to avoid visual ¡°zoom in/out¡± between frames.
%   ? Output format: GIF only. The frame interval equals FramePause.

% INPUT
%   Boundary      : struct for map boundary (same as mapMultipleCISLayout3D).
%                   Common fields: .X, .Y (vectors) for the study boundary.
%   System        : struct with the systems to draw (same schema as mapMultipleCISLayout3D).
%                   Typical fields:
%                     .PowerSystem / .GasSystem / .WaterSystem (each with .Node/.Edge arrays)
%                     .PowerToGas / .GasToPower / .PowerToWater / .WaterToPower (interdependency)
%   CascadeTrace  : 1¡ÁR struct array, each element describes one frame in time order.
%                   Required fields depend on Mode:
%                     Mode = 'Flow':
%                       (per system)  <Sys>EdgeState, <Sys>DemNode, <Sys>GenNode
%                       (links)        PowerToGas, GasToPower, PowerToWater, WaterToPower (optional)
%                     Mode = 'Connectivity':
%                       (per system)  <Sys>EdgeState, <Sys>NodeState
%                       (links)        PowerToGas, GasToPower, PowerToWater, WaterToPower (optional)
%                   Here <Sys> ¡Ê {'Power','Gas','Water'}.
%   Mode          : char, REQUIRED. Either 'Flow' or 'Connectivity'.
%                   ? 'Flow'        ¡ú calls mapMultipleCISLayout3D with Mode='Flow' and custom flows.
%                   ? 'Connectivity'¡ú calls mapMultipleCISLayout3D with Mode='Metric' and per-node/edge metrics.
%   options       : struct (minimal options consumed here; all others forwarded to the renderer)
%                   Consumed here:
%                     ? FigureColor : (default 'white') figure background color
%                     ? FramePause  : (default 1) seconds per frame; also used as GIF DelayTime
%                     ? GifFilename : (default 'CascadingProcess3D') base filename without extension;
%                                     if empty '', no file is written. The function saves '<GifFilename>.gif'.
%                   All other fields in `options` are passed through to mapMultipleCISLayout3D.
%                   For reliable animation, your mapMultipleCISLayout3D should respect:
%                     options.AxesOpt=hax (an axes handle)
% 
% OUTPUT
%   (none) This function renders an on-screen animation and optionally writes a GIF file.

%% Step 1: Required args & basic checks
if nargin < 4 || isempty(Mode)
    error('Mode is required and must be "Flow" or "Connectivity".');
end
if nargin < 5 || isempty(options), options = struct(); end

%% Step 2: Parse minimal animation options & prepare canvas
% Animation-level settings
[FigureColor, options] = popt(options,'FigureColor','white');
[FramePause,  options] = popt(options,'FramePause',1);
[GifFilename, options] = popt(options,'SaveGIF','CascadingProcess3D');

% Prepare figure/axes (created internally)
hfig = figure('Color', FigureColor, 'Name','Cascade Animation');
hax  = axes('Parent', hfig);
set(hax, 'NextPlot','replacechildren');

%% Step 3: Validate frames & initialize GIF
% Frames
R = numel(CascadeTrace);
if R==0
    warning('CascadeTrace is empty. Nothing to animate.');
    return;
end

% GIF init
doGIF = ~isempty(GifFilename);
delay = max(0.01, FramePause);   % GIF DelayTime

% For consistent view: lock after first frame
axesLocked = false;
v = struct();
sysList = {'Power','Gas','Water'};
ipsList = {'PowerToGas','GasToPower','PowerToWater','WaterToPower'};

%% Step 4: Main animation loop
for r = 1:R
    % Forward all remaining options to renderer:
    opts = options;  
    % build per-frame metrics from CascadeTrace
    switch Mode
        case 'Flow'
            PaintMode='Flow'; opts.FlowSource='custom';
            for s = 1:numel(sysList)
                if isfield(System,[sysList{s} 'System']) && ~isempty(System.([sysList{s} 'System']))
                    opts.([sysList{s} 'EdgeFlow']) = CascadeTrace(r).([sysList{s} 'EdgeState']);
                    opts.([sysList{s} 'NodeDem']) = CascadeTrace(r).([sysList{s} 'DemNode']);
                    opts.([sysList{s} 'NodeGen']) = CascadeTrace(r).([sysList{s} 'GenNode']);
                end
            end
            for i = 1:numel(ipsList)
                if isfield(System,ipsList{i}) && ~isempty(System.(ipsList{i})) && isfield(CascadeTrace(r),ipsList{i})
                    opts.([ipsList{i} 'Flow']) = CascadeTrace(r).(ipsList{i});
                end
            end
        case 'Connectivity'
            PaintMode='Metric'; 
            for s = 1:numel(sysList)
                if isfield(System,[sysList{s} 'System']) && ~isempty(System.([sysList{s} 'System']))
                    opts.([sysList{s} 'NodeMetric']) = CascadeTrace(r).([sysList{s} 'NodeState']);
                end
            end
            for i = 1:numel(ipsList)
                if isfield(System,ipsList{i}) && ~isempty(System.(ipsList{i})) && isfield(CascadeTrace(r),ipsList{i})
                    opts.([ipsList{i} 'Metric']) = CascadeTrace(r).(ipsList{i});
                end
            end
    end
    
    % Draw one frame
    cla(hax);
    opts.AxesOpt = hax;        % always draw on the same axes
    mapMultipleCISLayout3D(Boundary, System, PaintMode, opts);

    % Lock axis range after first frame so later frames won't rescale
    if ~axesLocked
        v.x = xlim(hax);  v.y = ylim(hax);  v.z = zlim(hax);
        v.da = daspect(hax);  v.pba = pbaspect(hax);
        v.cp = campos(hax);   v.ct = camtarget(hax);
        v.cva = camva(hax);   v.cu = camup(hax);
        v.proj = get(hax,'Projection');
        set(hax,'XLimMode','manual','YLimMode','manual','ZLimMode','manual', ...
            'DataAspectRatioMode','manual','PlotBoxAspectRatioMode','manual', ...
            'CameraViewAngleMode','manual', ...
            'XLim',v.x,'YLim',v.y,'ZLim',v.z);
        axesLocked = true;
    else
        set(hax,'XLim',v.x,'YLim',v.y,'ZLim',v.z);
        daspect(hax, v.da);  pbaspect(hax, v.pba);
        campos(hax, v.cp);   camtarget(hax, v.ct);  camva(hax, v.cva);  camup(hax, v.cu);
        set(hax,'Projection', v.proj);
    end

    drawnow;

    % Write GIF only
    if doGIF
        frame = getframe(hfig);
        [A,map] = rgb2ind(frame2im(frame),256);
        if r==1
            imwrite(A,map,[GifFilename '.gif'],'LoopCount',inf,'DelayTime',delay);
        else
            imwrite(A,map,[GifFilename '.gif'],'WriteMode','append','DelayTime',delay);
        end
    end
    pause(FramePause);
end
end 


function [val, s] = popt(s, name, default)
if isfield(s,name) && ~isempty(s.(name))
    val = s.(name);
    s = rmfield(s,name);
else
    val = default;
end
end