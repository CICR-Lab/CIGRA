function mapSingleCISRestorationDynamic(Boundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, options)
% INTRODUCTION
%   Visualize a single CIS's post-disaster restoration as a three-panel animation and export to GIF.
%   Panels:
%     (A) System functionality over time (left-continuous, piecewise-constant stairs)
%     (B) Network layout with damaged components highlighted; failures disappear as repairs finish
%     (C) Terminal-zone service levels across critical time points (choropleth by zone)
%
%   This version renders to a fixed-size figure (no output variables) and optionally writes a GIF.
%
% INPUTS
%   Boundary       : struct, optional map boundary for layout panel (B)
%                    - Fields (optional): X, Y (vectors or cell arrays for outer boundary)
%   CIS            : struct with CIS.Node(1¡ÁN), CIS.Edge(1¡ÁE)
%   TerminalZone   : 1¡ÁZ struct array describing terminal/service zones (fields: X, Y)
%   SysFunsEvo     : (Ndmg+1)¡Á4 numeric matrix: [time, normDrop, postFun, preFun]
%                    - Row 1 is pre-repair state at time = 0
%                    - Functionality ratio F(t) = postFun ./ preFun
%                    - Between t_i and t_{i+1}, F holds the value at t_i (stairs-left)
%   RepairSeq      : Ndmg¡Á5 numeric matrix: [DamageType, ComponentID, SysType, FinishTime, TeamID]
%                    - DamageType: 1=node, 2=edge; component is FAILED at time T iff FinishTime > T
%   ZoneStateEvo   : Z¡ÁK numeric matrix of zone states at critical time points
%                    - Column 1 is zone ID placeholder (not drawn); columns 2..K are service levels
%   options        : struct of optional parameters (all fields optional)
%     FigureColor       : default 'white'
%     CanvasWidth       : figure width in px (default 1280)
%     CanvasHeight      : figure height in px (default 360)
%     HideColorbars     : true/false, remove colorbars in (B) & (C) (default false)
%     EdgeWidth         : edge line width for layout (default 1.2)
%     PauseSec          : pause per frame for on-screen preview (default 1, only if ShowFigure=true)
%     TimeUnit          : time unit label for overlay text (default 'days')
%     TimeTextLocation  : {'northwest','northeast','southwest','southeast'} (default 'northeast')
%     TimeTextOnAxes    : true/false, draw time string on panel (A) axis (default true)
%     WriteGIF          : true/false, export GIF (default true)
%     GIFFile           : output filename (default 'CIS_Restoration.gif')
%     GIFDelaySec       : GIF frame delay in seconds (default 1)
%     ShowFigure        : true/false, preview on screen while rendering (default true)
%     % Layout color settings (optional):
%     MetricColormap    : colormap (M¡Á3) for layout panel (B) (default internal gradient)
%     ZoneColormap      : colormap (M¡Á3) for zone panel (C) (default internal gradient)
%     BoundaryColor     : RGB for boundary outline in (B) if used (default [0 0 0])
%     SystemType        : string passed to mapSingleCISLayout if needed (default [])
%
% OUTPUTS
%   (none) ¨C the function draws to a figure and (optionally) writes a GIF to disk.

%% Step 1: Validate inputs
if nargin < 7 || isempty(options), options = struct(); end

N = numel(CIS.Node); E = numel(CIS.Edge); Z = numel(TerminalZone);
if size(SysFunsEvo,2) < 4
    error('SysFunsEvo must be a (¡¤¡Á4) matrix: [time, normDrop, postFun, preFun].');
end
if ~isempty(RepairSeq) && size(RepairSeq,2) < 5
    error('RepairSeq must be (Ndmg¡Á5): [DamageType, ComponentID, SysType, FinishTime, TeamID].');
end
if size(ZoneStateEvo,1) ~= Z
    error('ZoneStateEvo row count (%d) must equal numel(TerminalZone)=%d.', size(ZoneStateEvo,1), Z);
end

%% Step 2: Parse options
% Options
FigureColor      = gopt(options,'FigureColor','white');
CanvasWidth      = gopt(options,'CanvasWidth',1280);
CanvasHeight     = gopt(options,'CanvasHeight',360);
HideColorbars    = gopt(options,'HideColorbars',true);

EdgeWidth        = gopt(options,'EdgeWidth',1.2);
PauseSec         = gopt(options,'PauseSec',1); % used only if ShowFigure=true
TimeUnit         = gopt(options,'TimeUnit','days');
TimeTextLocation = lower(gopt(options,'TimeTextLocation','northeast'));
TimeTextOnAxes   = gopt(options,'TimeTextOnAxes',true);
WriteGIF         = gopt(options,'WriteGIF',true);
GIFFile          = gopt(options,'GIFFile','CIS_Restoration.gif');
GIFDelaySec      = gopt(options,'GIFDelaySec',1);
ShowFigure       = gopt(options,'ShowFigure',true);

% Default colormaps
n = 256;
c = [0.18 0.23 1.00; 0.55 0.25 0.85; 1.00 0.34 0.20];
x = linspace(0,1,3);
xi = linspace(0,1,n);
cmapA = interp1(x,c,xi,'linear');
cmapA = flipud(cmapA*0.9 + 0.05);
cmapC=flipud([linspace(197,230,n)',linspace(224,80,n)',linspace(180,30,n)']./256);
defaultMetricMap = cmapA;
defaultZoneMap = cmapC;

%% Step 3: Build timeline & functionality
Times = SysFunsEvo(:,1).';
% Normalized functionality at each tiem 
funRatio = SysFunsEvo(:,3)./ SysFunsEvo(:,4); 
funRatio(~isfinite(funRatio)) = 0; 
funRatio = max(0,min(1,funRatio));

%% Step 4: Create figure & axes
fig = figure('Color',FigureColor, 'Name','Single CIS Restoration Dynamic', ...
    'Visible', tern(ShowFigure,'on','off'), 'Position',[100 100 CanvasWidth CanvasHeight]);

% Fixed axes layout to avoid jitter and reflow
axA = axes('Parent',fig,'Units','normalized','Position',[0.06 0.12 0.26 0.76], ...
    'ActivePositionProperty','position'); hold(axA,'on'); grid(axA,'on');
axB = axes('Parent',fig,'Units','normalized','Position',[0.38 0.12 0.28 0.76], ...
    'ActivePositionProperty','position'); hold(axB,'on');
axC = axes('Parent',fig,'Units','normalized','Position',[0.70 0.12 0.26 0.76], ...
    'ActivePositionProperty','position'); hold(axC,'on');

xlabel(axA,'Time'); ylabel(axA,'System functionality');

% Fixed annotation titles (annotation avoids axis-title jitter)
place_title_annotation(fig, axA, '(A) Function over time');
place_title_annotation(fig, axB, '(B) Network damage & repair');
place_title_annotation(fig, axC, '(C) Zone service');

% Time overlay placement (axes-based so every frame is captured)
switch TimeTextLocation
    case 'northwest',  txpos = [0.02, 0.96]; ha='left';  va='top';
    case 'northeast',  txpos = [0.98, 0.96]; ha='right'; va='top';
    case 'southwest',  txpos = [0.02, 0.04]; ha='left';  va='bottom';
    case 'southeast',  txpos = [0.98, 0.04]; ha='right'; va='bottom';
    otherwise,         txpos = [0.98, 0.96]; ha='right'; va='top';
end
if TimeTextOnAxes
    txh = text(axA, txpos(1), txpos(2), '', 'Units','normalized', ...
        'HorizontalAlignment',ha, 'VerticalAlignment',va, 'FontWeight','bold', ...
        'Interpreter','none', 'Clipping','off');
else
    % fallback figure-level annotation
    txh = annotation(fig,'textbox','String','', 'Units','normalized', ...
        'Position',[txpos(1) txpos(2) 0.12 0.06], 'FitBoxToText','on', ...
        'HorizontalAlignment',ha, 'VerticalAlignment',va, 'EdgeColor','none', ...
        'BackgroundColor','none','FontWeight','bold');
end

%% Step 5: Preplot reference & initialize
[xs_full, ys_full] = stairs_like(Times, funRatio);
plot(axA, xs_full, ys_full, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
lnAnim = plot(axA, xs_full(1), ys_full(1), '-', 'LineWidth', 2.0);
xlim(axA, [min(Times) max(Times)] + 0.02*[-1 1]*max(1,range(Times)) );
ylim(axA, [0 1]);

%% Step 6: Animation loop
for it = 1:numel(Times)
    T = Times(it);
    
% (A) update stairs up to current time T
    [xA,yA] = stairs_upto(Times, funRatio, T);
    set(lnAnim,'XData',xA,'YData',yA);
    
% update time overlay text (recreate if deleted)
    if ~exist('txh','var') || ~isgraphics(txh)
        if TimeTextOnAxes
            txh = text(axA, txpos(1), txpos(2), '', 'Units','normalized', ...
                'HorizontalAlignment',ha, 'VerticalAlignment',va, 'FontWeight','bold', ...
                'Interpreter','none', 'Clipping','off');
        else
            txh = annotation(fig,'textbox','String','', 'Units','normalized', ...
                'Position',[txpos(1) txpos(2) 0.12 0.06], 'FitBoxToText','on', ...
                'HorizontalAlignment',ha, 'VerticalAlignment',va, 'EdgeColor','none', ...
                'BackgroundColor','none','FontWeight','bold');
        end
    end
    if isa(txh,'matlab.graphics.primitive.Text') || isa(txh,'matlab.graphics.shape.TextBox')
        txh.String = sprintf('%s %s', sprintf('t = %.0f', T), TimeUnit);
    end
    
% (B) draw layout with damage status at time T
    cla(axB); hold(axB,'on');
    [NodeMetric, EdgeMetric] = damage_metric_at_time(N, E, RepairSeq, T);
    optB = struct();
    optB.SystemType    = gopt(options,'SystemType',[]);
    optB.EdgeMetric    = EdgeMetric(:);
    optB.NodeMetric    = NodeMetric(:);
    optB.FlowColormap  = gopt(options,'MetricColormap',defaultMetricMap);
    optB.MetricLimits  = [0,1];
    optB.EdgeWidth     = EdgeWidth;
    optB.BoundaryColor = gopt(options,'BoundaryColor',[0 0 0]);
    optB.FigureColor   = FigureColor;
    optB.AxesOpt       = axB;
    optB.NodeMinSize   = 20;
    
    mapSingleCISLayout(Boundary, CIS, 'Metric', optB);
    if HideColorbars
        delete(findall(fig,'Type','ColorBar'));
    end
    
% (C) draw zone snapshot at time T
    cla(axC); hold(axC,'on');
    zVals = zone_values_at_time(ZoneStateEvo, Times, T);
    optC = struct('AxesOpt', axC, 'Colormap', gopt(options,'ZoneColormap',defaultZoneMap), ...
        'FigureColor', FigureColor,'ValueLimits',[0,1]);
    mapTerminalZoneValue(TerminalZone, zVals, optC);
    if HideColorbars
        delete(findall(fig,'Type','ColorBar'));
    end
    
    drawnow;
    
% write GIF frame if requested
    if WriteGIF
        fr = getframe(fig);
        [im, map] = frame2im(fr);
        if isempty(map)
            [imind, map] = rgb2ind(im, 256, 'nodither');
        else
            imind = rgb2ind(im, map);
        end
        if it == 1
            imwrite(imind, map, GIFFile, 'gif', 'LoopCount', inf, 'DelayTime', GIFDelaySec);
        else
            imwrite(imind, map, GIFFile, 'gif', 'WriteMode', 'append', 'DelayTime', GIFDelaySec);
        end
    end
    
    if ShowFigure && PauseSec > 0
        pause(PauseSec);
    end
end
end

function th = place_title_annotation(fig, ax, str)
pos = ax.Position; % normalized
padY = 0.02; % small gap above the axes
th = annotation(fig,'textbox','String',str, 'Units','normalized', ...
    'Position',[pos(1) pos(2)+pos(4)+padY pos(3) 0.04], ...
    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
    'EdgeColor','none','BackgroundColor','none','FontWeight','bold');
end

function out = tern(cond,a,b)
if cond, out=a; else, out=b; end
end

function [x,y] = stairs_like(t, v)
% Construct left-continuous stairs that holds value at t_i on [t_i, t_{i+1})
[t, idx] = unique(t(:),'stable'); v = v(idx);
if numel(t) == 1, x=[t;t]; y=[v;v]; return; end
x = [t(1); reshape([t(2:end)'; t(2:end)'], [], 1)];
y = [v(1); reshape([v(1:end-1)'; v(2:end)'], [], 1)];
x = [t(1); x]; y = [v(1); y];
end

function [x,y] = stairs_upto(t, v, T)
% Stairs truncated at time T (left-continuous)
[t, idx] = unique(t(:),'stable'); v = v(idx);
mask = t <= T; lastIdx = find(mask,1,'last');
if isempty(lastIdx)
    x = t(1); y = v(1); return;
end
% Use values up to lastIdx, then extend horizontal segment to T
x = t(1); y = v(1);
for i = 2:lastIdx
    x = [x; t(i); t(i)]; %#ok<AGROW>
    y = [y; y(end); v(i)]; %#ok<AGROW>
end
x = [x; T]; y = [y; y(end)];
end

function [NodeMetric, EdgeMetric] = damage_metric_at_time(N, E, RepairSeq, T)
% NodeMetric/EdgeMetric are 1 for FAILED, 0 for HEALTHY at time T
NodeMetric = ones(N,1); EdgeMetric = ones(E,1);
if isempty(RepairSeq), return; end
DT  = RepairSeq(:,1);  % 1=node, 2=edge
CID = RepairSeq(:,2);
FT  = RepairSeq(:,4);

% Failed iff FinishTime > T (i.e., not yet repaired)
notRepaired = FT > T;
DT = DT(notRepaired); CID = CID(notRepaired);

% Clamp to valid ranges
isNode = DT == 1; isEdge = DT == 2;
nodeIDs = CID(isNode);
edgeIDs = CID(isEdge); 

NodeMetric(nodeIDs) = 0;
EdgeMetric(edgeIDs) = 0;
end


function zVals = zone_values_at_time(ZoneStateEvo, funTime, T)
% Returns a Z¡Á1 vector of zone values at time T by selecting the latest
% critical time point <= T. If ZoneTimes is empty, align to funTime index.
Z = size(ZoneStateEvo,1);
K = size(ZoneStateEvo,2); % includes ID column

% Align to functionality timeline length (excluding time=0 safety)
t = funTime(:)';
idx = find(t <= T, 1, 'last');
if isempty(idx), idx = 1; end
% Cap to available zone columns (2..K)
idx = max(1, min(idx, K-1));
zVals = ZoneStateEvo(:, idx+1); % shift by 1 due to ID column
end

function v = gopt(s,name,default)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function out = range(x)
if isempty(x), out = 0; else, out = max(x)-min(x); end
end
