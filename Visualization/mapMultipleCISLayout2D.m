function mapMultipleCISLayout2D(Boundary, System, Mode, options)
% INTRODUCTION:
%   mapMultipleCISLayout2D renders up to three CIS layers (Power/Gas/Water)
%   on a 2D map by shifting each layer along one axis ('x' or 'y').
%   It supports three visualization modes:
%     - 'Topology' - structural view: edges (polylines) and nodes (points).
%     - 'Flow'     - per-edge flow (color) + per-node generation/demand (size/type);
%                    data can come from system fields or user-provided custom vectors.
%     - 'Metric'   - per-edge / per-node user metrics (e.g., 0/1 states, failure probability).
%
% INPUTS:
%   Boundary : struct with optional fields:
%              - X, Y : equal-length coordinates of one continuous study
%                boundary. A trailing NaN is allowed. Used to draw a colored
%                stripe per layer.
%
%   System   : struct that may include any of the following systems (any subset allowed):
%              - PowerSystem.Node(1-by-Np), .Edge(1-by-Ep)
%              - GasSystem.Node(1-by-Ng),   .Edge(1-by-Eg)
%              - WaterSystem.Node(1-by-Nw), .Edge(1-by-Ew)
%              Node must have at least: Longitude, Latitude.
%              Edge is recommended to have: X, Y (NaN-separated polyline).
%              Optional interdependency tables (rows = links):
%              - PowerToGas, GasToPower, PowerToWater, WaterToPower
%                (assumed columns include [fromNode, toNode, ...]; if custom Flow/Metric is provided
%                 via options, it will be used to color links in non-Topology modes.)
%
%   Mode     : 'Topology' | 'Flow' | 'Metric'
%
%   options  : struct controlling styles, data sources and layout:
%     General (all modes)
%       - FigureColor      : default 'white'
%       - BoundaryColor    : default [0 0 0]
%       - PlaneFaceAlpha   : default 0.30           % stripe fill transparency
%       - NetworkColor     : default [96 96 96]/255 % neutral gray
%       - MarkerSizeTopo   : default 4
%       - EdgeWidth        : default 1
%       - ThemeColorPower  : default [255 187 0]/255 % amber
%       - ThemeColorGas    : default [124 187 0]/255 % green
%       - ThemeColorWater  : default [1 159 212]/255 % cyan
%       - LayoutAxis       : 'x' or 'y' (default 'y') % choose offset direction
%       - LayerGap         : inter-layer gap; default = 5% of span along LayoutAxis
%       - AxesOpt          : handle to existing axes for drawing
%       - ValueLabel       : colorbar label (default 'Normalized |flow|' in
%                            Flow mode and 'Metric' in Metric mode)
%
%     Flow mode
%       - FlowSource       : 'system' (default) | 'custom'
%       - FlowColormap     : M-by-3 numeric in [0,1]; default parula(256)
%       - NodeMinSize/NodeMaxSize : node size range (default 10/20; also used in Metric mode)
%       - SourceFaceColor  : color for source/generator nodes (default [0.15 0.47 0.80])
%       - DemandFaceColor  : color for demand nodes (default [0.86 0.33 0.10])
%       - IpsLinkColor     : M-by-3 colormap for interdependency links
%       - When FlowSource='custom', the following are REQUIRED for every plotted system:
%           PowerEdgeFlow (Ep-by-1), PowerNodeGen (Np-by-1), PowerNodeDem (Np-by-1)
%           GasEdgeFlow   (Eg-by-1), GasNodeGen   (Ng-by-1), GasNodeDem   (Ng-by-1)
%           WaterEdgeFlow (Ew-by-1), WaterNodeGen (Nw-by-1), WaterNodeDem (Nw-by-1)
%       - Optional custom interdependency flow magnitude for coloring links:
%           PowerToGasFlow, GasToPowerFlow, PowerToWaterFlow, WaterToPowerFlow
%         (vector length must equal the corresponding table's row count)
%
%     Metric mode
%       - FlowColormap     : M-by-3 numeric in [0,1]; default = a soft red-blue map (256)
%       - MetricLimits     : optional [min max] shared color limits
%       - Per-system metric vectors are optional:
%           PowerEdgeMetric (Ep-by-1), PowerNodeMetric (Np-by-1)
%           GasEdgeMetric   (Eg-by-1), GasNodeMetric   (Ng-by-1)
%           WaterEdgeMetric (Ew-by-1), WaterNodeMetric (Nw-by-1)
%         Missing metrics use NetworkColor. Provided metrics share one color scale.
%       - Optional custom interdependency metric for coloring links:
%           PowerToGasMetric, GasToPowerMetric, PowerToWaterMetric, WaterToPowerMetric
%
% OUTPUT:
%   The function draws on options.AxesOpt when supplied; otherwise it creates
%   a new figure and axes:
%     - A 2D map with up to three horizontally/vertically offset stripes (by LayoutAxis).
%     - In 'Flow' and 'Metric' modes, a shared colorbar is added.
%     - Interdependency links are drawn between corresponding node pairs; in non-Topology
%       modes they are colored using provided Flow/Metric (normalized to [0,1]) if available.
%   No variables are returned.

if nargin < 3 || isempty(Mode)
    error('Argument "mode" is required and must be "Topology", "Flow" or "Metric".');
end
if nargin < 4 || isempty(options), options = struct(); end

Mode = string(Mode);
validModes = ["Topology","Flow","Metric"];
if ~any(Mode == validModes)
    error('Invalid mode: %s. Use "Topology", "Flow" or "Metric".', Mode);
end

FigureColor     = gopt(options,'FigureColor','white');
BoundaryColor   = gopt(options,'BoundaryColor',[0 0 0]);
PlaneFaceAlpha  = gopt(options,'PlaneFaceAlpha',0.30);
NetworkColor    = gopt(options,'NetworkColor',[96 96 96]/255);
MarkerSizeTopo  = gopt(options,'MarkerSizeTopo',4);
EdgeWidth       = gopt(options,'EdgeWidth',1);

FlowColormap    = gopt(options,'FlowColormap',[]);
NodeMinSize     = gopt(options,'NodeMinSize',10);
NodeMaxSize     = gopt(options,'NodeMaxSize',20);
SourceFaceColor = gopt(options,'SourceFaceColor',[0.15 0.47 0.80]);
DemandFaceColor = gopt(options,'DemandFaceColor',[0.86 0.33 0.10]);
IpsLinkColor    = gopt(options,'IpsLinkColor',[]);
ThemeColorPower = gopt(options,'ThemeColorPower',[255 187 0]/255); % amber
ThemeColorGas   = gopt(options,'ThemeColorGas',[124 187 0]/255); % green
ThemeColorWater = gopt(options,'ThemeColorWater',[1 159 212]/255); % cyan
LayoutAxis      = gopt(options,'LayoutAxis','y');  % 'x' or 'y'
AxesOpt         = gopt(options,'AxesOpt', []);
ValueLabel      = gopt(options,'ValueLabel',[]);

CLim            = gopt(options,'MetricLimits', []);

LayoutAxis = lower(string(LayoutAxis));
if ~isscalar(LayoutAxis) || ~any(LayoutAxis == ["x","y"])
    error('options.LayoutAxis must be "x" or "y".');
end

[xmin,xmax,ymin,ymax] = localBounds(Boundary, System);
if strcmp(LayoutAxis,'y')
    stripeSpan = ymax - ymin;
else
    stripeSpan = xmax - xmin;
end
if ~isfinite(stripeSpan) || stripeSpan <= 0
    stripeSpan = 1;
end
defaultGap = 0.05 * stripeSpan;
gap = gopt(options,'LayerGap', defaultGap);
if ~isfinite(gap) || gap < 0, gap = defaultGap; end
LayerStep = stripeSpan + gap;

nColors = 256;
c = [0.18 0.23 1.00; 0.55 0.25 0.85; 1.00 0.34 0.20];
x = linspace(0,1,3);
xi = linspace(0,1,nColors);
cmapA = interp1(x,c,xi,'linear');
cmapA = flipud(cmapA*0.9 + 0.05);
cmapC = parula(256);
if strcmp(Mode,'Metric')
    defaultFlowMap = cmapA;
else
    defaultFlowMap = cmapC;
end

flowmap = resolveColormap(FlowColormap, 'FlowColormap', defaultFlowMap);
ipsmap = resolveColormap(IpsLinkColor, 'IpsLinkColor', cmapA);

hasPower = isfield(System,'PowerSystem') && ~isempty(System.PowerSystem);
hasGas   = isfield(System,'GasSystem')   && ~isempty(System.GasSystem);
hasWater = isfield(System,'WaterSystem') && ~isempty(System.WaterSystem);
hasSystem=[hasPower,hasGas,hasWater];
sysList = {'Power','Gas','Water'};
linkFields = {'PowerToGas','GasToPower','PowerToWater','WaterToPower'};
Interdependency = struct();

if ~any(hasSystem)
    warning('No systems found to draw.');
    return;
end

ordP = 2;
ordG = 2 + hasPower;
ordW = 1 + (~hasPower && ~hasGas);
[offP,offG,offW] = deal(NaN);
if hasPower, offP = (ordP-1)*LayerStep; end
if hasGas,   offG = (ordG-1)*LayerStep; end
if hasWater, offW = (ordW-1)*LayerStep; end

FlowSource = string(gopt(options,'FlowSource','system'));
if strcmp(Mode,'Flow')
    if ~any(FlowSource == ["system","custom"])
        error('options.FlowSource must be "system" or "custom" when mode="Flow".');
    end
end

linkSuffix = '';
if strcmp(Mode,'Flow') && strcmp(FlowSource,'custom')
    for s = 1:numel(sysList)
        if hasSystem(s)
            N = numel(System.([sysList{s} 'System']).Node); E = numel(System.([sysList{s} 'System']).Edge);
            EdgeFlow.(sysList{s}) = getOptVec(options,[sysList{s} 'EdgeFlow'],E, 0, [sysList{s} 'EdgeFlow'],true);
            NodeGen.(sysList{s})  = getOptVec(options,[sysList{s} 'NodeGen'], N, 0, [sysList{s} 'NodeGen'],true);
            NodeDem.(sysList{s})  = getOptVec(options,[sysList{s} 'NodeDem'], N, 0, [sysList{s} 'NodeDem'],true);
        end
    end
    linkSuffix = 'Flow';
end

if strcmp(Mode,'Metric')
    for s = 1:numel(sysList)
        if hasSystem(s)
            N = numel(System.([sysList{s} 'System']).Node); E = numel(System.([sysList{s} 'System']).Edge);
            EdgeMetric.(sysList{s}) = getOptVec(options,[sysList{s} 'EdgeMetric'],E, NaN, [sysList{s} 'EdgeMetric'],false);
            NodeMetric.(sysList{s}) = getOptVec(options,[sysList{s} 'NodeMetric'],N, NaN, [sysList{s} 'NodeMetric'],false);
        end
    end
    linkSuffix = 'Metric';
    if ~isempty(CLim)
        if ~isnumeric(CLim) || numel(CLim) ~= 2 || ~isreal(CLim) || any(~isfinite(CLim(:)))
            error('options.MetricLimits must be a finite 1-by-2 numeric vector [min max].');
        end
        CLim = CLim(:).';
        if CLim(1) >= CLim(2)
            error('options.MetricLimits must satisfy MetricLimits(1) < MetricLimits(2).');
        end
    end
    metricVals = [];
    for s = 1:numel(sysList)
        if hasSystem(s)
            em = EdgeMetric.(sysList{s}); nm = NodeMetric.(sysList{s});
            metricVals = [metricVals; em(:); nm(:)]; %#ok<AGROW>
        end
    end
    if isempty(CLim), metricLims = pickLims(metricVals); else, metricLims = CLim; end
end

if ~isempty(linkSuffix)
    for linkIdx = 1:numel(linkFields)
        field = linkFields{linkIdx};
        optionName = [field linkSuffix];
        if isfield(System,field) && ~isempty(System.(field)) && ...
                isfield(options,optionName) && ~isempty(options.(optionName))
            Interdependency.(field) = [System.(field), ...
                getOptVec(options,optionName,size(System.(field),1),NaN,optionName,false)];
        end
    end
end

if ~isempty(AxesOpt)
    ax  = AxesOpt;
else
    fig = figure('Color', FigureColor);
    ax  = axes('Parent', fig);
end
holdState = ishold(ax);
set(ax,'color','none'); hold(ax,'on'); grid(ax,'on');
axis(ax,'equal');

if isfield(Boundary,'X') && isfield(Boundary,'Y')
    if hasPower, drawBoundaryStripe(Boundary, offP, ThemeColorPower); end
    if hasGas,   drawBoundaryStripe(Boundary, offG, ThemeColorGas);   end
    if hasWater, drawBoundaryStripe(Boundary, offW, ThemeColorWater); end
end

switch Mode
    case 'Topology'
        drawTopologyLayer('power');
        drawTopologyLayer('gas');
        drawTopologyLayer('water');
    case 'Flow'
        drawFlowLayer('power');
        drawFlowLayer('gas'  );
        drawFlowLayer('water');
        colormap(ax, flowmap); caxis(ax,[0 1]); cb = colorbar('peer',ax);
        if isempty(ValueLabel), ValueLabel = 'Normalized |flow|'; end
        cb.Label.String = ValueLabel; cb.Label.FontSize = 14; cb.FontSize = 14;
    case 'Metric'
        drawMetricLayer('power');
        drawMetricLayer('gas'  );
        drawMetricLayer('water');
        colormap(ax, flowmap); caxis(ax,metricLims); cb = colorbar('peer',ax);
        if isempty(ValueLabel), ValueLabel = 'Metric'; end
        cb.Label.String = ValueLabel; cb.Label.FontSize = 14; cb.FontSize = 14;
end

plotLinks('PowerToGas','power',offP,'gas',offG);
plotLinks('GasToPower','gas',offG,'power',offP);
plotLinks('PowerToWater','power',offP,'water',offW);
plotLinks('WaterToPower','water',offW,'power',offP);

axis(ax,'off');
if ~holdState, hold(ax,'off'); end

    function drawBoundaryStripe(B, off, faceColor)
        lon = B.X(:)'; lat = B.Y(:)';
        [x,y] = applyOffset(lon,lat,off);
        plot(ax,x,y,'-','Color',BoundaryColor,'LineWidth',EdgeWidth);
        m = ~isnan(x) & ~isnan(y);
        if any(m)
            fill(ax,x(m), y(m), faceColor, 'FaceAlpha', PlaneFaceAlpha, 'EdgeColor','none');
        end
    end

    function drawTopologyLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                off = offP; col = ThemeColorPower; S = System.PowerSystem;
            case 'gas'
                if ~hasGas, return; end
                off = offG; col = ThemeColorGas; S = System.GasSystem;
            case 'water'
                if ~hasWater, return; end
                off = offW; col = ThemeColorWater; S = System.WaterSystem;
        end
        
        E = S.Edge;
        for e = 1:numel(E)
            [ex,ey] = getXY(E(e)); if isempty(ex), continue; end
            [ex,ey] = applyOffset(ex,ey,off);
            plot(ax,ex,ey,'-','Color',col,'LineWidth',EdgeWidth);
        end
        
        N = S.Node;
        for nodeIdx = 1:numel(N)
            [x,y] = applyOffset(N(nodeIdx).Longitude,N(nodeIdx).Latitude,off);
            td = fget(N(nodeIdx),'TargetDemand',0);
            if td > 0, mkr = 'o'; else, mkr = 's'; end
            plot(ax,x,y,mkr,'MarkerSize',MarkerSizeTopo,'Color',col,'MarkerFaceColor',col);
        end
    end

    function drawFlowLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                off = offP; S = System.PowerSystem;
            case 'gas'
                if ~hasGas, return; end
                off = offG; S = System.GasSystem;
            case 'water'
                if ~hasWater, return; end
                off = offW; S = System.WaterSystem;
        end
        if strcmp(FlowSource,'custom')
            switch kind
                case 'power', flowVals = EdgeFlow.Power; genVals = NodeGen.Power; demVals = NodeDem.Power;
                case 'gas', flowVals = EdgeFlow.Gas; genVals = NodeGen.Gas; demVals = NodeDem.Gas;
                case 'water', flowVals = EdgeFlow.Water; genVals = NodeGen.Water; demVals = NodeDem.Water;
            end
        else
            flowVals = arrayfun(@(e) fget(e,'RealFlow',NaN), S.Edge);
            genVals = arrayfun(@(node) fget(node,'RealGeneration',0), S.Node);
            demVals = arrayfun(@(node) fget(node,'RealDemand',0), S.Node);
        end
        flowForColor = abs(flowVals);
        flowForColor(~isfinite(flowForColor)) = NaN;
        genLims  = pickLims(genVals);
        demLims  = pickLims(demVals);
        flowMaxField = [upper(kind(1)) kind(2:end) 'FlowMax'];
        flowScale = gopt(options,flowMaxField,[]);
        if ~isempty(flowScale) && (~isnumeric(flowScale) || ~isscalar(flowScale) || ~isfinite(flowScale) || flowScale <= 0)
            error('options.%s must be a positive finite scalar.',flowMaxField);
        end
        if isempty(flowScale)
            finiteFlow = flowForColor(isfinite(flowForColor));
            if isempty(finiteFlow), flowScale = 1; else, flowScale = max(finiteFlow); end
        end
        if flowScale <= 0, flowScale = 1; end
        flowNorm = flowForColor / flowScale;
        E = S.Edge;
        for e = 1:numel(E)
            [ex,ey] = getXY(E(e)); if isempty(ex), continue; end
            [ex,ey] = applyOffset(ex,ey,off);
            c = value2color(flowNorm(e), [0 1], flowmap, NetworkColor);
            plot(ax,ex,ey,'-','Color',c,'LineWidth',EdgeWidth);
        end
        
        N = S.Node;
        isSrc = genVals>0;
        isDem = (~isSrc) & (demVals>0);
        isOth = (~isSrc) & (~isDem);
        gSize = mapToSize(genVals, genLims, NodeMinSize, NodeMaxSize);
        dSize = mapToSize(demVals, demLims, NodeMinSize, NodeMaxSize);
        oSize = ones(size(genVals))*NodeMinSize*0.8;
        
        if any(isSrc)
            [x,y] = nxys(N,isSrc); [x,y]=applyOffset(x,y,off);
            scatter(ax,x,y,gSize(isSrc),'s','filled','MarkerFaceColor',SourceFaceColor,'MarkerFaceAlpha',0.95);
        end
        if any(isDem)
            [x,y] = nxys(N,isDem); [x,y]=applyOffset(x,y,off);
            scatter(ax,x,y,dSize(isDem),'o','filled','MarkerFaceColor',DemandFaceColor,'MarkerFaceAlpha',0.95);
        end
        if any(isOth)
            [x,y] = nxys(N,isOth); [x,y]=applyOffset(x,y,off);
            scatter(ax,x,y,oSize(isOth),'^','filled','MarkerFaceColor',NetworkColor,'MarkerFaceAlpha',0.95);
        end
    end

    function drawMetricLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                off = offP;  S=System.PowerSystem;
            case 'gas'
                if ~hasGas, return; end
                off = offG;  S=System.GasSystem;
            case 'water'
                if ~hasWater, return; end
                off = offW;  S=System.WaterSystem;
        end
        switch kind
            case 'power', EM = EdgeMetric.Power; NM = NodeMetric.Power;
            case 'gas', EM = EdgeMetric.Gas; NM = NodeMetric.Gas;
            case 'water', EM = EdgeMetric.Water; NM = NodeMetric.Water;
        end
        for e = 1:numel(S.Edge)
            [ex,ey] = getXY(S.Edge(e)); if isempty(ex), continue; end
            col = value2color(EM(e), metricLims, flowmap, NetworkColor);
            [ex,ey] = applyOffset(ex,ey,off);
            plot(ax,ex,ey,'-','Color',col,'LineWidth',EdgeWidth);
        end
        [x_all,y_all] = nxys(S.Node, true(1,numel(S.Node)));
        [x,y] = applyOffset(x_all,y_all,off);
        C = zeros(numel(S.Node),3);
        for i=1:numel(S.Node)
            C(i,:) = value2color(NM(i), metricLims, flowmap, NetworkColor);
        end
        sSize = mapToSize(NM, metricLims, NodeMinSize, NodeMaxSize);
        scatter(ax, x, y, sSize, C, 'o', 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.95);
    end

    function plotLinks(field, fromK, offFrom, toK, offTo)
        if ~isfield(System, field) || isempty(System.(field)), return; end
        L = System.(field);
        if ~strcmp(Mode,'Topology')
            if strcmp(Mode,'Flow') && strcmp(FlowSource,'system')
                ipsVals = ones(size(L,1),1);
            elseif isfield(Interdependency,field)
                linkData = Interdependency.(field);
                if strcmp(Mode,'Flow')
                    if strcmp(field,'GasToPower'), denomCol = 5; else, denomCol = 3; end
                    ipsVals = linkData(:,end)./linkData(:,denomCol);
                else
                    ipsVals = linkData(:,end);
                end
            else
                ipsVals = NaN(size(L,1),1);
            end
        end
        for i = 1:size(L,1)
            [fx,fy] = nodeCoord2D(fromK, L(i,1), offFrom);
            [tx,ty] = nodeCoord2D(toK,   L(i,2), offTo);
            if ~strcmp(Mode,'Topology')
                col = value2color(ipsVals(i), [0,1], ipsmap, NetworkColor);
            else
                col =NetworkColor;
            end
            plot(ax,[fx tx],[fy ty],'-','Color',col,'LineWidth',EdgeWidth);
        end
    end

    function [x,y] = applyOffset(x,y,off)
        if strcmp(LayoutAxis,'x')
            x = x + off; % shift horizontally
        else
            y = y + off; % shift vertically
        end
    end

    function [x,y] = nodeCoord2D(kind, idx, off)
        switch kind
            case 'power'
                N = System.PowerSystem.Node;
            case 'gas'
                N = System.GasSystem.Node;
            case 'water'
                N = System.WaterSystem.Node;
        end
        if idx>0 && idx<=numel(N)
            x = N(idx).Longitude; y = N(idx).Latitude; [x,y] = applyOffset(x,y,off);
        else
            x = NaN; y = NaN;
        end
    end
end


function [xmin,xmax,ymin,ymax] = localBounds(Boundary, System)
xs = []; ys = [];
if isfield(Boundary,'X'), xs = [xs; Boundary.X(:)]; ys = [ys; Boundary.Y(:)]; end
if isfield(System,'PowerSystem') && isfield(System.PowerSystem,'Node')
    xs = [xs; [System.PowerSystem.Node.Longitude]'];
    ys = [ys; [System.PowerSystem.Node.Latitude]'];
end
if isfield(System,'GasSystem') && isfield(System.GasSystem,'Node')
    xs = [xs; [System.GasSystem.Node.Longitude]'];
    ys = [ys; [System.GasSystem.Node.Latitude]'];
end
if isfield(System,'WaterSystem') && isfield(System.WaterSystem,'Node')
    xs = [xs; [System.WaterSystem.Node.Longitude]'];
    ys = [ys; [System.WaterSystem.Node.Latitude]'];
end
xs = xs(isfinite(xs)); ys = ys(isfinite(ys));
if isempty(xs) || isempty(ys)
    xmin=0; xmax=1; ymin=0; ymax=1;
else
    xmin=min(xs); xmax=max(xs); ymin=min(ys); ymax=max(ys);
end
end

function M = resolveColormap(C, name, fallback)
if isempty(C)
    M = fallback; return;
end
if ~isnumeric(C) || ~ismatrix(C) || size(C,2)~=3 || ~isreal(C)
    error('Invalid %s: must be a nonempty M-by-3 numeric matrix with values in [0,1].', name);
end
if any(~isfinite(C(:))) || any(C(:)<0 | C(:)>1)
    error('%s must contain finite values in [0,1].', name);
end
if size(C,1) == 1
    M = repmat(C,256,1);
elseif size(C,1) < 256
    M = resampleColormap(C,256);
else
    M = C;
end
end

function cmap2 = resampleColormap(cmap, targetN)
m = size(cmap,1);
cmap2 = interp1(1:m,cmap,linspace(1,m,targetN),'linear');
end


function v = gopt(s, name, default)
if isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function [ex,ey] = getXY(edge)
ex = []; ey = [];
if isfield(edge,'X') && isfield(edge,'Y')
    ex = edge.X; ey = edge.Y;
end
end

function val = fget(s, field, default)
if isfield(s, field) && ~isempty(s.(field)), val = full(s.(field)); else, val = default; end
end

function c = value2color(v, lims, cmap, fallback)
if ~isfinite(v)
    c = fallback; return;
end
v = max(min(v, lims(2)), lims(1));
t = (v - lims(1)) / (lims(2) - lims(1));
k = 1 + floor(t*(size(cmap,1)-1));
c = cmap(k,:);
end

function sz = mapToSize(vals, lims, smin, smax)
vals(~isfinite(vals)) = NaN;
vals = abs(vals);
sz = (vals - lims(1)) / (lims(2)-lims(1));
sz = max(0, min(1, sz));
sz = smin + sz .* (smax - smin);
sz(~isfinite(sz)) = smin;
end

function lims = pickLims(v)
w = v(isfinite(v));
if isempty(w)
    lims = [0 1];
else
    lo = min(w); hi = max(w);
    if hi <= lo, hi = lo + eps; end
    lims = [lo hi];
end
end

function [x,y] = nxys(nodes,mask)
nodes = nodes(mask);
x = [nodes.Longitude].';
y = [nodes.Latitude].';
end

function v = getOptVec(s, name, n, defaultVal, label, required)
if n == 0
    v = zeros(0,1);
    return;
end
if isfield(s,name) && ~isempty(s.(name))
    v0 = s.(name)(:);
    if ~isnumeric(v0) || numel(v0) ~= n
        error('%s length (%d) does not match expected length (%d).', label, numel(v0), n);
    else
        v = v0;
    end
else
    if required
        error('options.%s is required and must contain %d values.',label,n);
    end
    v = repmat(defaultVal,n,1);
end
end
