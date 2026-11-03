function mapMultipleCISLayout3D(Boundary, System, Mode, options)
% INTRODUCTION:
%   This function renders multiple CIS layers (Power/Gas/Water) on separate z-planes and
%   optionally draws interdependency links between layers. It supports three modes:
%     - 'Topology' - Structural view (nodes + edges) per system.
%     - 'Flow'     - Edge flow magnitudes (color) + node generation/demand (size/type).
%                     Flows can come from system fields or user-provided custom inputs.
%     - 'Metric'   - User-defined scalar metrics for each system (per-edge/per-node),
%                     e.g., failure probability (0-1), binary states, importance scores.
%
% INPUTS:
%   Boundary  : struct with optional X and Y coordinates for one continuous
%               study boundary. A trailing NaN is allowed.
%   System    : struct with up to three systems and optional interdependency tables:
%               - System.PowerSystem.Node(1-by-Np), .Edge(1-by-Ep)
%               - System.GasSystem.Node(1-by-Ng),   .Edge(1-by-Eg)
%               - System.WaterSystem.Node(1-by-Nw), .Edge(1-by-Ew)
%                 Each Node needs at least Longitude, Latitude. Each Edge ideally has
%                 X, Y polyline coordinates (NaN-separated) for proper geometry.
%               - Optional interdependency tables (examples):
%                 System.PowerToGas, GasToPower, PowerToWater, WaterToPower
%                 (assumed columns: [fromNode, toNode, cap, flow, ...]).
%
%   Mode      : 'Topology' | 'Flow' | 'Metric'
%
%   options   : struct to control styles, data sources, and metrics.
%               General (all modes):
%                 - FigureColor      ('white')
%                 - BoundaryColor    ([0 0 0])
%                 - PlaneFaceAlpha   (0.30)      % filled plane transparency
%                 - NetworkColor     ([96 96 96]/255)
%                 - MarkerSizeTopo   (4)
%                 - EdgeWidth        (1)
%                 - ThemeColorPower  ([255 187 0]/255)   % amber
%                 - ThemeColorGas    ([124 187 0]/255)   % green
%                 - ThemeColorWater  ([1 159 212]/255)   % cyan
%                 - AxesOpt          : handle to existing axes for drawing
%                 - ValueLabel       : colorbar label (default 'Normalized |flow|' in
%                                      Flow mode and 'Metric' in Metric mode)
%
%               Flow mode:
%                 - FlowSource       :'system' (default) | 'custom'
%                 - FlowColormap     : M-by-3 numeric in [0,1] (default parula(256))
%                 - IpsLinkColor     : M-by-3 interdependency-link colormap
%                 - NodeMinSize/NodeMaxSize : node size range (also used in Metric mode)
%                 - SourceFaceColor  : [0.15 0.47 0.80]
%                 - DemandFaceColor  : [0.86 0.33 0.10]
%                 - Custom inputs (required when FlowSource='custom'):
%                      PowerEdgeFlow (Ep-by-1), PowerNodeGen (Np-by-1), PowerNodeDem (Np-by-1)
%                      GasEdgeFlow   (Eg-by-1), GasNodeGen   (Ng-by-1), GasNodeDem   (Ng-by-1)
%                      WaterEdgeFlow (Ew-by-1), WaterNodeGen (Nw-by-1), WaterNodeDem (Nw-by-1)
%                 - Optional custom interdependency (when plotting links with values):
%                      PowerToGasFlow, GasToPowerFlow, PowerToWaterFlow, WaterToPowerFlow
%
%               Metric mode:
%                 - FlowColormap     : M-by-3 numeric in [0,1]; default = a soft red-blue map (256)
%                 - MetricLimits     : optional [min max] shared color limits
%                 - Edge/Node metrics (optional; missing values use NetworkColor):
%                      PowerEdgeMetric (Ep-by-1), PowerNodeMetric (Np-by-1)
%                      GasEdgeMetric   (Eg-by-1), GasNodeMetric   (Ng-by-1)
%                      WaterEdgeMetric (Ew-by-1), WaterNodeMetric (Nw-by-1)
%                 - Optional custom interdependency metric (for link colors):
%                      PowerToGasMetric, GasToPowerMetric, PowerToWaterMetric, WaterToPowerMetric
%
% OUTPUT:
%   The function draws on options.AxesOpt when supplied; otherwise it creates
%   a new figure and axes:
%     - A 3D map with up to three systems on separate z-planes.
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
AxesOpt         = gopt(options,'AxesOpt', []);
ValueLabel      = gopt(options,'ValueLabel',[]);

CLim            = gopt(options,'MetricLimits', []);

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
z_power = 2;
z_gas = 2 + hasPower;
z_water = 1 + (~hasPower && ~hasGas);

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
set(ax,'color','none'); hold(ax,'on'); grid(ax,'on'); view(ax,3);
axis(ax,'equal');

if isfield(Boundary,'X') && isfield(Boundary,'Y')
    lon = Boundary.X(:)'; lat = Boundary.Y(:)';
    m = ~isnan(lon) & ~isnan(lat);
    if hasPower, drawBoundaryPlane(z_power,ThemeColorPower); end
    if hasGas, drawBoundaryPlane(z_gas,ThemeColorGas); end
    if hasWater, drawBoundaryPlane(z_water,ThemeColorWater); end
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
        colormap(ax,flowmap); caxis(ax,[0 1]); cb = colorbar('peer',ax);
        if isempty(ValueLabel), ValueLabel = 'Normalized |flow|'; end
        cb.Label.String = ValueLabel; cb.Label.FontSize = 14; cb.FontSize = 14;
    case 'Metric'
        drawMetricLayer('power');
        drawMetricLayer('gas'  );
        drawMetricLayer('water');
        colormap(ax,flowmap); caxis(ax,metricLims); cb = colorbar('peer',ax);
        if isempty(ValueLabel), ValueLabel = 'Metric'; end
        cb.Label.String = ValueLabel; cb.Label.FontSize = 14; cb.FontSize = 14;
end

plotLinks('PowerToGas',  'power','gas');
plotLinks('GasToPower',  'gas',  'power');
plotLinks('PowerToWater','power','water');
plotLinks('WaterToPower','water','power');

axis(ax,'off'); view(ax,63,33);
if ~holdState, hold(ax,'off'); end

    function drawBoundaryPlane(z, faceColor)
        plot3(ax,lon,lat,z*ones(size(lon)),'-','Color',BoundaryColor,'LineWidth',EdgeWidth);
        fill3(ax,lon(m),lat(m),z*ones(1,sum(m)),faceColor, ...
            'FaceAlpha',PlaneFaceAlpha,'EdgeColor','none');
    end

    function drawTopologyLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                z = z_power; col = ThemeColorPower; S = System.PowerSystem;
            case 'gas'
                if ~hasGas, return; end
                z = z_gas; col = ThemeColorGas; S = System.GasSystem;
            case 'water'
                if ~hasWater, return; end
                z = z_water; col = ThemeColorWater; S = System.WaterSystem;
        end
        
        E = S.Edge;
        for e = 1:numel(E)
            [ex,ey] = getXY(E(e)); if isempty(ex), continue; end
            plot3(ax,ex,ey,z*ones(size(ex)),'-','Color',col,'LineWidth',EdgeWidth);
        end
        
        N = S.Node;
        for nodeIdx = 1:numel(N)
            td = fget(N(nodeIdx),'TargetDemand',0);
            if td > 0, mkr = 'o'; else, mkr = 's'; end
            plot3(ax,N(nodeIdx).Longitude,N(nodeIdx).Latitude,z,mkr,'MarkerSize',MarkerSizeTopo, ...
                'Color',col,'MarkerFaceColor',col);
        end
    end

    function drawFlowLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                z = z_power; S = System.PowerSystem;
            case 'gas'
                if ~hasGas, return; end
                z = z_gas; S = System.GasSystem;
            case 'water'
                if ~hasWater, return; end
                z = z_water; S = System.WaterSystem;
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
            c = value2color(flowNorm(e), [0 1], flowmap, NetworkColor);
            plot3(ax,ex,ey,z*ones(size(ex)),'-','Color',c,'LineWidth',EdgeWidth);
        end
        
        N = S.Node;
        isSrc = genVals>0;
        isDem = (~isSrc) & (demVals>0);
        isOth = (~isSrc) & (~isDem);
        gSize = mapToSize(genVals, genLims, NodeMinSize, NodeMaxSize);
        dSize = mapToSize(demVals, demLims, NodeMinSize, NodeMaxSize);
        oSize = ones(size(genVals))*NodeMinSize*0.8;
        
        if any(isSrc)
            [x,y] = nxys(N,isSrc);
            scatter3(ax,x,y,z*ones(size(x)), gSize(isSrc), 's', 'filled', ...
                'MarkerFaceColor', SourceFaceColor,  'MarkerFaceAlpha', 0.95);
        end
        if any(isDem)
            [x,y] = nxys(N,isDem);
            scatter3(ax,x,y,z*ones(size(x)), dSize(isDem), 'o', 'filled', ...
                'MarkerFaceColor', DemandFaceColor, 'MarkerFaceAlpha', 0.95);
        end
        if any(isOth)
            [x,y] = nxys(N,isOth);
            scatter3(ax,x,y,z*ones(size(x)), oSize(isOth), '^', 'filled', ...
                'MarkerFaceColor', NetworkColor, 'MarkerFaceAlpha', 0.95);
        end
    end

    function drawMetricLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                z = z_power; S=System.PowerSystem;
            case 'gas'
                if ~hasGas, return; end
                z = z_gas; S=System.GasSystem;
            case 'water'
                if ~hasWater, return; end
                z = z_water; S=System.WaterSystem;
        end
        switch kind
            case 'power', EM = EdgeMetric.Power; NM = NodeMetric.Power;
            case 'gas', EM = EdgeMetric.Gas; NM = NodeMetric.Gas;
            case 'water', EM = EdgeMetric.Water; NM = NodeMetric.Water;
        end
        for e = 1:numel(S.Edge)
            [ex,ey] = getXY(S.Edge(e)); if isempty(ex), continue; end
            col = value2color(EM(e), metricLims, flowmap, NetworkColor);
            plot3(ax,ex,ey,z*ones(size(ex)),'-','Color',col,'LineWidth',EdgeWidth);
        end
        [x_all,y_all] = nxys(S.Node, true(1,numel(S.Node)));
        C = zeros(numel(S.Node),3);
        for i=1:numel(S.Node)
            C(i,:) = value2color(NM(i), metricLims, flowmap, NetworkColor);
        end
        sSize = mapToSize(NM, metricLims, NodeMinSize, NodeMaxSize);
        scatter3(ax, x_all, y_all, z*ones(1,numel(S.Node)), sSize, C, 'o', 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.95);
    end

    function plotLinks(field, fromK, toK)
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
            [fx,fy,fz] = nodeCoord(fromK, L(i,1));
            [tx,ty,tz] = nodeCoord(toK,   L(i,2));
            if ~strcmp(Mode,'Topology')
                col = value2color(ipsVals(i), [0,1], ipsmap, NetworkColor);
            else
                col =NetworkColor;
            end
            plot3(ax,[fx tx],[fy ty],[fz tz],'-','Color',col,'LineWidth',EdgeWidth);
        end
    end

    function [x,y,z] = nodeCoord(kind, idx)
        switch kind
            case 'power'
                N = System.PowerSystem.Node; z = z_power;
            case 'gas'
                N = System.GasSystem.Node;  z = z_gas;
            case 'water'
                N = System.WaterSystem.Node; z = z_water;
        end
        if idx>0 && idx<=numel(N)
            x = N(idx).Longitude; y = N(idx).Latitude;
        else
            x = NaN; y = NaN; z = NaN;
        end
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
    v = repmat(defaultVal, n, 1);
end
end
