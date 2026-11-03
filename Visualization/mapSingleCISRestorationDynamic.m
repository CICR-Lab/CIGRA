function mapSingleCISRestorationDynamic(Boundary, CIS, TerminalZone, SysFunsEvo, RepairSeq, ZoneStateEvo, options)
% INTRODUCTION
%   Visualize a single CIS's post-disaster restoration as a three-panel
%   animation and optionally export it to GIF.
%   Panels:
%     (A) System functionality over time (piecewise-constant stairs)
%     (B) Network layout with damaged components highlighted; failures
%         disappear as repairs finish
%     (C) Terminal-zone service levels at critical time points
%
% INPUTS
%   Boundary       : struct, optional map boundary for layout panel (B)
%                    - Fields (optional): X, Y (vectors or cell arrays for outer boundary)
%   CIS            : struct with CIS.Node(1-by-N), CIS.Edge(1-by-E)
%   TerminalZone   : 1-by-Z struct array describing terminal/service zones
%                    (fields: X, Y)
%   SysFunsEvo     : (Ndmg+1)-by-4 numeric matrix:
%                    [time, normDrop, postFun, preFun]
%                    - Row 1 is pre-repair state at time = 0
%                    - Functionality ratio F(t) = postFun ./ preFun
%                    - Between t_i and t_{i+1}, F holds the value at t_i
%   RepairSeq      : Ndmg-by-5 numeric matrix:
%                    [DamageType, ComponentID, SysType, FinishTime, TeamID]
%                    - DamageType: 1=node, 2=edge; a component is failed at
%                      time T when FinishTime > T
%   ZoneStateEvo   : Z-by-(size(SysFunsEvo,1)+1) numeric matrix
%                    - Column 1 contains zone IDs; remaining columns contain
%                      service levels aligned with rows of SysFunsEvo
%   options        : struct containing required SystemType and the optional
%                    parameters below
%     SystemType        : 'power' | 'gas' | 'water' | 'road' (required)
%     FigureColor       : default 'white'
%     CanvasWidth       : figure width in px (default 1280)
%     CanvasHeight      : figure height in px (default 360)
%     HideColorbars     : true/false, remove colorbars in (B) & (C) (default true)
%     EdgeWidth         : edge line width for layout (default 1.2)
%     PauseSec          : pause per frame for on-screen preview (default 0.2,
%                         only if ShowFigure=true)
%     TimeUnit          : time unit label for overlay text (default 'days')
%     TimeTextLocation  : {'northwest','northeast','southwest','southeast'} (default 'northeast')
%     TimeTextOnAxes    : true/false, draw time text on panel (A) axes or as
%                         a figure annotation beside panel (A) (default true)
%     WriteGIF          : true/false, export GIF (default true)
%     GIFFile           : output filename (default 'CIS_Restoration.gif')
%     GIFDelaySec       : GIF frame delay in seconds (default 1)
%     ShowFigure        : true/false, preview on screen while rendering (default true)
%     MetricColormap    : M-by-3 colormap for layout panel (B)
%                         (default internal gradient)
%     ZoneColormap      : M-by-3 colormap for zone panel (C)
%                         (default internal gradient)
%     BoundaryColor     : RGB for boundary outline in (B) if used (default [0 0 0])
%
% OUTPUTS
%   None. The function draws to a figure and optionally writes a GIF.

%% Step 1: Validate inputs
if nargin < 7 || isempty(options), options = struct(); end

SystemType = string(gopt(options,'SystemType',''));
if ~isscalar(SystemType) || strlength(SystemType) == 0
    error('options.SystemType is required: power | gas | water | road.');
end
if ~any(SystemType == ["power","gas","water","road"])
    error('options.SystemType must be one of: power, gas, water, road.');
end

N = numel(CIS.Node); E = numel(CIS.Edge); Z = numel(TerminalZone);
if size(SysFunsEvo,2) < 4
    error('SysFunsEvo must have at least 4 columns: [time, normDrop, postFun, preFun].');
end
if ~isempty(RepairSeq) && size(RepairSeq,2) < 5
    error('RepairSeq must have at least 5 columns: [DamageType, ComponentID, SysType, FinishTime, TeamID].');
end
if size(ZoneStateEvo,1) ~= Z
    error('ZoneStateEvo row count (%d) must equal numel(TerminalZone)=%d.', size(ZoneStateEvo,1), Z);
end
expectedZoneCols = size(SysFunsEvo,1) + 1;
if size(ZoneStateEvo,2) ~= expectedZoneCols
    error('ZoneStateEvo must have %d columns: one zone-ID column plus one state column per SysFunsEvo row.', expectedZoneCols);
end

%% Step 2: Parse options
FigureColor      = gopt(options,'FigureColor','white');
CanvasWidth      = gopt(options,'CanvasWidth',1280);
CanvasHeight     = gopt(options,'CanvasHeight',360);
HideColorbars    = gopt(options,'HideColorbars',true);

EdgeWidth        = gopt(options,'EdgeWidth',1.2);
PauseSec         = gopt(options,'PauseSec',0.2); % used only if ShowFigure=true
TimeUnit         = gopt(options,'TimeUnit','days');
TimeTextLocation = lower(gopt(options,'TimeTextLocation','northeast'));
TimeTextOnAxes   = gopt(options,'TimeTextOnAxes',true);
WriteGIF         = gopt(options,'WriteGIF',true);
GIFFile          = gopt(options,'GIFFile','CIS_Restoration.gif');
GIFDelaySec      = gopt(options,'GIFDelaySec',1);
ShowFigure       = gopt(options,'ShowFigure',true);

n = 256;
c = [0.18 0.23 1.00; 0.55 0.25 0.85; 1.00 0.34 0.20];
x = linspace(0,1,3);
xi = linspace(0,1,n);
cmapA = interp1(x,c,xi,'linear');
cmapA = flipud(cmapA*0.9 + 0.05);
cmapC=flipud([linspace(197,230,n)',linspace(224,80,n)',linspace(180,30,n)']./256);
defaultMetricMap = cmapA;
defaultZoneMap = cmapC;
MetricColormap = gopt(options,'MetricColormap',defaultMetricMap);
ZoneColormap   = gopt(options,'ZoneColormap',defaultZoneMap);
BoundaryColor  = gopt(options,'BoundaryColor',[0 0 0]);

%% Step 3: Build timeline & functionality
StateTimes = SysFunsEvo(:,1).';
funRatio = SysFunsEvo(:,3)./ SysFunsEvo(:,4); 
funRatio(~isfinite(funRatio)) = 0; 
funRatio = max(0,min(1,funRatio));
% At simultaneous completion times, keep the final state and draw one frame.
[~, lastIdx] = unique(StateTimes,'last');
lastIdx = sort(lastIdx);
Times = StateTimes(lastIdx);
funRatio = funRatio(lastIdx);

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

xlabel(axA,sprintf('Time (%s)',TimeUnit)); ylabel(axA,'System functionality');

% Fixed annotation titles (annotation avoids axis-title jitter)
place_title_annotation(fig, axA, '(A) Function over time');
place_title_annotation(fig, axB, '(B) Network damage & repair');
place_title_annotation(fig, axC, '(C) Zone service');

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
    annPos = annotation_position(axA, txpos, ha, va);
    txh = annotation(fig,'textbox','String','', 'Units','normalized', ...
        'Position',annPos, 'FitBoxToText','off', ...
        'HorizontalAlignment',ha, 'VerticalAlignment',va, 'EdgeColor','none', ...
        'BackgroundColor','none','FontWeight','bold');
end

%% Step 5: Preplot reference & initialize
[xs_full, ys_full] = stairs(Times, funRatio);
plot(axA, xs_full, ys_full, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
lnAnim = plot(axA, xs_full(1), ys_full(1), '-', 'LineWidth', 2.0);
xlim(axA, [min(Times) max(Times)] + 0.02*[-1 1]*max(1,max(Times)-min(Times)) );
ylim(axA, [0 1]);

%% Step 6: Animation loop
for it = 1:numel(Times)
    T = Times(it);

    [xA,yA] = stairs_upto(Times, funRatio, T);
    set(lnAnim,'XData',xA,'YData',yA);
    txh.String = sprintf('t = %g %s', T, TimeUnit);

    cla(axB);
    [NodeMetric, EdgeMetric] = damage_metric_at_time(N, E, RepairSeq, T);
    optB = struct('SystemType',SystemType,'EdgeMetric',EdgeMetric, ...
        'NodeMetric',NodeMetric,'FlowColormap',MetricColormap, ...
        'MetricLimits',[0,1],'ValueLabel','Component functionality', ...
        'EdgeWidth',EdgeWidth,'BoundaryColor',BoundaryColor, ...
        'AxesOpt',axB,'NodeMinSize',20);
    mapSingleCISLayout(Boundary, CIS, 'Metric', optB);

    cla(axC);
    zVals = zone_values_at_time(ZoneStateEvo, StateTimes, T);
    optC = struct('AxesOpt', axC, 'Colormap', ZoneColormap, ...
        'ValueLimits',[0,1], 'ValueLabel','Zone functionality');
    mapTerminalZoneValue(TerminalZone, zVals, optC);
    if HideColorbars
        delete(findall(fig,'Type','ColorBar'));
    end

    drawnow;

    if WriteGIF
        im = frame2im(getframe(fig));
        [imind, map] = rgb2ind(im, 256, 'nodither');
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

function place_title_annotation(fig, ax, str)
pos = ax.Position;
padY = 0.02;
annotation(fig,'textbox','String',str, 'Units','normalized', ...
    'Position',[pos(1) pos(2)+pos(4)+padY pos(3) 0.04], ...
    'HorizontalAlignment','center','VerticalAlignment','bottom', ...
    'EdgeColor','none','BackgroundColor','none','FontWeight','bold');
end

function pos = annotation_position(ax, txpos, ha, va)
pos = [ax.Position(1) + txpos(1)*ax.Position(3), ...
    ax.Position(2) + txpos(2)*ax.Position(4), 0.12, 0.06];
if strcmp(ha,'right'), pos(1) = pos(1) - pos(3); end
if strcmp(va,'top'), pos(2) = pos(2) - pos(4); end
end

function out = tern(cond,a,b)
if cond, out=a; else, out=b; end
end

function [x,y] = stairs_upto(t, v, T)
t = t(:); v = v(:);
lastIdx = find(t <= T,1,'last');
[x,y] = stairs(t(1:lastIdx),v(1:lastIdx));
if x(end) < T
    x(end+1) = T;
    y(end+1) = y(end);
end
end

function [NodeMetric, EdgeMetric] = damage_metric_at_time(N, E, RepairSeq, T)
% NodeMetric/EdgeMetric are 0 for failed and 1 for operational at time T.
NodeMetric = ones(N,1); EdgeMetric = ones(E,1);
if isempty(RepairSeq), return; end
DT  = RepairSeq(:,1);  % 1=node, 2=edge
CID = RepairSeq(:,2);
FT  = RepairSeq(:,4);

% Failed iff FinishTime > T (i.e., not yet repaired)
notRepaired = FT > T;
DT = DT(notRepaired); CID = CID(notRepaired);

isNode = DT == 1; isEdge = DT == 2;
nodeIDs = CID(isNode);
edgeIDs = CID(isEdge); 

NodeMetric(nodeIDs) = 0;
EdgeMetric(edgeIDs) = 0;
end


function zVals = zone_values_at_time(ZoneStateEvo, funTime, T)
idx = find(funTime <= T,1,'last');
zVals = ZoneStateEvo(:,idx+1);
end

function v = gopt(s,name,default)
if isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end
