function mapMultipleCISLayout2D(Boundary, System, Mode, options)
% INTRODUCTION:
%   mapMultipleCISLayout2D renders up to three CIS layers (Power/Gas/Water)
%   on a 2D map by shifting each layer along one axis ('x' or 'y').
%   It supports three visualization modes:
%     - 'Topology' ¡ª structural view: edges (polylines) and nodes (points).
%     - 'Flow'     ¡ª per-edge flow (color) + per-node generation/demand (size/type);
%                    data can come from system fields or user-provided custom vectors.
%     - 'Metric'   ¡ª per-edge / per-node user metrics (e.g., 0/1 states, failure probability).
%
% INPUTS:
%   Boundary : struct with optional fields:
%              - X, Y : study boundary (polyline; NaN breaks allowed). Used to draw a colored stripe per layer.
%
%   System   : struct that may include any of the following systems (any subset allowed):
%              - PowerSystem.Node(1¡ÁNp), .Edge(1¡ÁEp)
%              - GasSystem.Node(1¡ÁNg),   .Edge(1¡ÁEg)
%              - WaterSystem.Node(1¡ÁNw), .Edge(1¡ÁEw)
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
%
%     Flow mode
%       - FlowSource       : 'system' (default) | 'custom'
%       - FlowColormap     : M¡Á3 numeric in [0,1]; default parula(256)
%       - NodeMinSize/NodeMaxSize : node size range (default 10/20)
%       - SourceFaceColor  : color for source/generator nodes (default [0.15 0.47 0.80])
%       - DemandFaceColor  : color for demand nodes (default [0.86 0.33 0.10])
%       - When FlowSource='custom', the following are REQUIRED for every plotted system:
%           PowerEdgeFlow (Ep¡Á1), PowerNodeGen (Np¡Á1), PowerNodeDem (Np¡Á1)
%           GasEdgeFlow   (Eg¡Á1), GasNodeGen   (Ng¡Á1), GasNodeDem   (Ng¡Á1)
%           WaterEdgeFlow (Ew¡Á1), WaterNodeGen (Nw¡Á1), WaterNodeDem (Nw¡Á1)
%       - Optional custom interdependency magnitude for coloring links:
%           PowerToGasFlow, GasToPowerFlow, PowerToWaterFlow, WaterToPowerFlow
%         (vector length must equal the corresponding table¡¯s row count)
%
%     Metric mode
%       - FlowColormap     : M¡Á3 numeric in [0,1]; default = a soft red¨Cblue map (256)
%       - For each plotted system, the following are REQUIRED:
%           PowerEdgeMetric (Ep¡Á1), PowerNodeMetric (Np¡Á1)
%           GasEdgeMetric   (Eg¡Á1), GasNodeMetric   (Ng¡Á1)
%           WaterEdgeMetric (Ew¡Á1), WaterNodeMetric (Nw¡Á1)
%         (values can be in [0,1] for probabilities, or any scalar domain; a shared colorbar is used)
%       - Optional custom interdependency metric for coloring links:
%           PowerToGasMetric, GasToPowerMetric, PowerToWaterMetric, WaterToPowerMetric
%
% OUTPUT:
%   The function draws directly on the current axes:
%     - A 2D map with up to three horizontally/vertically offset stripes (by LayoutAxis).
%     - In 'Flow' and 'Metric' modes, a shared colorbar is added.
%     - Interdependency links are drawn between corresponding node pairs; in non-Topology
%       modes they are colored using provided Flow/Metric (normalized to [0,1]) if available.
%   No variables are returned.

%% Step 1: required args & basic checks
if nargin < 3 || isempty(Mode)
    error('Argument "mode" is required and must be "Topology", "Flow" or "Metric".');
end
if nargin < 4 || isempty(options), options = struct(); end

% Validate mode value
Mode = string(Mode);
validModes = ["Topology","Flow","Metric"];
if ~any(Mode == validModes)
    error('Invalid mode: %s. Use "Topology", "Flow" or "Metric".', Mode);
end

%% Step 2: Validate inputs
% Load style / behavior options
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

CLim            = gopt(options,'MetricLimits', []);

% Compute layer spacing (stripe size + gap) from bounds
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

% Default colormaps
n = 256;
c = [0.18 0.23 1.00; 0.55 0.25 0.85; 1.00 0.34 0.20];
x = linspace(0,1,3);
xi = linspace(0,1,n);
cmapA = interp1(x,c,xi,'linear');
cmapA = flipud(cmapA*0.9 + 0.05);
cmapB = [linspace(0.90,0.40,n)' linspace(0.45,0.65,n)' linspace(0.45,0.95,n)'];
cmapC = parula(256);
if strcmp(Mode,'Metric')
    defaultFlowMap = cmapA;
else
    defaultFlowMap = cmapC;
end

% Resolve working colormaps
flowmap = resolveColormap(FlowColormap, 'FlowColormap', defaultFlowMap);
ipsmap = resolveColormap(IpsLinkColor, 'IpsLinkColor', cmapA);

% which systems exist
hasPower = isfield(System,'PowerSystem') && ~isempty(System.PowerSystem);
hasGas   = isfield(System,'GasSystem')   && ~isempty(System.GasSystem);
hasWater = isfield(System,'WaterSystem') && ~isempty(System.WaterSystem);
hasSystem=[hasPower,hasGas,hasWater];
sysList = {'Power','Gas','Water'};

% Layer order assignments
if hasPower && hasGas && hasWater
    ordP = 2; ordG = 3; ordW = 1;
elseif hasPower && hasGas
    ordP = 2; ordG = 3; ordW = NaN;
elseif hasPower && hasWater
    ordP = 2; ordG = NaN; ordW = 1;
elseif hasGas && hasWater
    ordP = NaN; ordG = 2; ordW = 1;
elseif hasPower
    ordP = 2;  ordG = NaN; ordW = NaN;
elseif hasGas
    ordP = NaN; ordG = 2;  ordW = NaN;
elseif hasWater
    ordP = NaN; ordG = NaN; ordW = 2;
end

% Convert order indices to actual offsets
[offP,offG,offW] = deal(NaN);
if hasPower, offP = (ordP-1)*LayerStep; end
if hasGas,   offG = (ordG-1)*LayerStep; end
if hasWater, offW = (ordW-1)*LayerStep; end

% FlowSource handling: only required/validated when mode=="Flow"
FlowSource = gopt(options,'FlowSource',[]);
if strcmp(Mode,'Flow')
    if isempty(FlowSource)
        warning('mode="Flow" without options.FlowSource. Defaulting to "system" (system-provided values).');
        FlowSource = "system";
    else
        FlowSource = string(FlowSource);
        if ~any(FlowSource == ["system","custom"])
            error('options.FlowSource must be "system" or "custom" when mode="Flow".');
        end
    end
else
    FlowSource = NaN;
end

% Validate custom inputs if requested
if strcmp(FlowSource,'custom') && strcmp(Mode,'Flow')
    for s = 1:numel(sysList)
        if hasSystem(s)==1
            N = numel(System.([sysList{s} 'System']).Node); E = numel(System.([sysList{s} 'System']).Edge);
            EdgeFlow.(sysList{s}) = getOptVec(options,[sysList{s} 'EdgeFlow'],E, 0, [sysList{s} 'EdgeFlow']);
            NodeGen.(sysList{s})  = getOptVec(options,[sysList{s} 'NodeGen'], N, 0, [sysList{s} 'NodeGen']);
            NodeDem.(sysList{s})  = getOptVec(options,[sysList{s} 'NodeDem'], N, 0, [sysList{s} 'NodeDem']);
        else
            EdgeFlow.(sysList{s}) = []; NodeGen.(sysList{s}) = []; NodeDem.(sysList{s}) = [];
        end
    end
    if hasPower&&hasGas
        Interdependency.PowerToGas=System.PowerToGas; Interdependency.GasToPower=System.GasToPower;
        Interdependency.PowerToGas(:,end+1)=getOptVec(options,'PowerToGasFlow',size(System.PowerToGas,1), 1, 'PowerToGasFlow');
        Interdependency.GasToPower(:,end+1)=getOptVec(options,'GasToPowerFlow',size(System.GasToPower,1), 1, 'GasToPowerFlow');
    end
    if hasPower&&hasWater
        Interdependency.PowerToWater=System.PowerToWater; Interdependency.WaterToPower=System.WaterToPower;
        Interdependency.PowerToWater(:,end+1)=getOptVec(options,'PowerToWaterFlow',size(System.PowerToWater,1), 1, 'PowerToWaterFlow');
        Interdependency.WaterToPower(:,end+1)=getOptVec(options,'WaterToPowerFlow',size(System.WaterToPower,1), 1, 'WaterToPowerFlow');
    end
end

% Require metrics when Mode='Metric'
if strcmp(Mode,'Metric')
    for s = 1:numel(sysList)
        if hasSystem(s)==1
            N = numel(System.([sysList{s} 'System']).Node); E = numel(System.([sysList{s} 'System']).Edge);
            EdgeMetric.(sysList{s}) = getOptVec(options,[sysList{s} 'EdgeMetric'],E, 1, [sysList{s} 'EdgeMetric']);
            NodeMetric.(sysList{s}) = getOptVec(options,[sysList{s} 'NodeMetric'],N, 0, [sysList{s} 'NodeMetric']);
        end
    end
    if hasPower&&hasGas
        Interdependency.PowerToGas=System.PowerToGas; Interdependency.GasToPower=System.GasToPower;
        Interdependency.PowerToGas(:,end+1)=getOptVec(options,'PowerToGasMetric',size(System.PowerToGas,1), 0, 'PowerToGasMetric');
        Interdependency.GasToPower(:,end+1)=getOptVec(options,'GasToPowerMetric',size(System.GasToPower,1), 0, 'GasToPowerMetric');
    end
    if hasPower&&hasWater
        Interdependency.PowerToWater=System.PowerToWater; Interdependency.WaterToPower=System.WaterToPower;
        Interdependency.PowerToWater(:,end+1)=getOptVec(options,'PowerToWaterMetric',size(System.PowerToWater,1), 0, 'PowerToWaterMetric');
        Interdependency.WaterToPower(:,end+1)=getOptVec(options,'WaterToPowerMetric',size(System.WaterToPower,1), 0, 'WaterToPowerMetric');
    end
    if ~isempty(CLim)
        if ~isnumeric(CLim) || numel(CLim) ~= 2
            error('options.MetricLimits must be a 1¡Á2 numeric vector [min max].');
        end
        if CLim(1) >= CLim(2)
            error('options.MetricLimits must satisfy MetricLimits(1) < MetricLimits(2).');
        end
    end
end

%% Step 3: Figure & Axes
if ~isempty(AxesOpt)
    ax  = options.AxesOpt;
    fig = ancestor(ax,'figure');
    if exist('FigureColor','var') && ~isempty(FigureColor)
        set(fig,'Color', FigureColor);
    end
else
    fig = figure('Color', FigureColor);
    ax  = axes('Parent', fig);
end
ax = gca; set(ax,'color','none'); hold(ax,'on'); grid(ax,'on');
axis equal;

% Draw boundary stripes per existing layer
if isfield(Boundary,'X') && isfield(Boundary,'Y')
    if hasPower, drawBoundaryStripe(Boundary, offP, ThemeColorPower); end
    if hasGas,   drawBoundaryStripe(Boundary, offG, ThemeColorGas);   end
    if hasWater, drawBoundaryStripe(Boundary, offW, ThemeColorWater); end
end

%% Step 4: Plot system by mode
% Dispatch to per-layer drawers
switch Mode
    case 'Topology'
        drawTopologyLayer('power');
        drawTopologyLayer('gas');
        drawTopologyLayer('water');
    case 'Flow'
        drawFlowLayer('power');
        drawFlowLayer('gas'  );
        drawFlowLayer('water');
        colormap(ax, flowmap); cb = colorbar('peer',ax); cb.Label.String ='|Flow|';
    case 'Metric'
        drawMetricLayer('power');
        drawMetricLayer('gas'  );
        drawMetricLayer('water');
        colormap(ax, flowmap); cb = colorbar('peer',ax); cb.Label.String ='Metric';
end

%% Step 5: Plot interdependent links
% Draw cross-layer links; in non-Topology modes color by normalized value.
plotLinks('PowerToGas','power',offP,'gas',offG);
plotLinks('GasToPower','gas',offG,'power',offP);
plotLinks('PowerToWater','power',offP,'water',offW);
plotLinks('WaterToPower','water',offW,'power',offP);

axis off; hold off;

%% helper functions
    function drawBoundaryStripe(B, off, faceColor)
        lon = B.X(:)'; lat = B.Y(:)';
        [x,y] = applyOffset(lon,lat,off);
        plot(x,y,'-','Color',BoundaryColor,'LineWidth',EdgeWidth);
        m = ~isnan(x) & ~isnan(y);
        if any(m)
            fill(x(m), y(m), faceColor, 'FaceAlpha', PlaneFaceAlpha, 'EdgeColor','none');
        end
    end

    function drawTopologyLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                off = offP; col = ThemeColorPower; S = System.PowerSystem; %%
            case 'gas'
                if ~hasGas, return; end
                off = offG; col = ThemeColorGas; S = System.GasSystem;
            case 'water'
                if ~hasWater, return; end
                off = offW; col = ThemeColorWater; S = System.WaterSystem;
        end
        
        % --- Edges ---
        E = S.Edge;
        for e = 1:numel(E)
            [ex,ey] = getXY(E(e)); if isempty(ex), continue; end
            [ex,ey] = applyOffset(ex,ey,off);
            plot(ex,ey,'-','Color',col,'LineWidth',EdgeWidth);
        end
        
        % --- Nodes ---
        N = S.Node;
        for n = 1:numel(N)
            [x,y] = applyOffset(N(n).Longitude,N(n).Latitude,off);
            td = fget(N(n),'TargetDemand',0);
            mkr = iff(td>0,'o','s');
            plot(x,y,mkr,'MarkerSize',MarkerSizeTopo,'Color',col,'MarkerFaceColor',col);
        end
    end

    function drawFlowLayer(kind)
        switch kind
            case 'power'
                if ~hasPower, return; end
                off = offP; S = System.PowerSystem;%%
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
            genVals = arrayfun(@(n) fget(n,'MaxGeneration',0), S.Node);
            demVals = arrayfun(@(n) fget(n,'TargetDemand',0), S.Node);
        end
        % --- Normalize for color ---
        flowForColor = abs(flowVals);
        flowForColor(~isfinite(flowForColor)) = NaN;
        % --- Limits  ---
        flowLims = pickLims(flowForColor);
        genLims  = pickLims(genVals);
        demLims  = pickLims(demVals);
        % --- edges colored by flow ---
        E = S.Edge;
        for e = 1:numel(E)
            [ex,ey] = getXY(E(e)); if isempty(ex), continue; end
            [ex,ey] = applyOffset(ex,ey,off);
            c = value2color(flowForColor(e), flowLims, flowmap, NetworkColor);
            plot(ex,ey,'-','Color',c,'LineWidth',EdgeWidth);
        end
        if all(isfinite(flowLims)) && diff(flowLims) > 0, caxis(ax, flowLims); end
        
        % --- nodes sized & colored by type ---
        N = S.Node;
        isSrc = genVals>0;
        isDem = (~isSrc) & (demVals>0);
        isOth = (~isSrc) & (~isDem);
        gSize = mapToSize(genVals, genLims, NodeMinSize, NodeMaxSize);
        dSize = mapToSize(demVals, demLims, NodeMinSize, NodeMaxSize);
        oSize = ones(size(genVals))*NodeMinSize*0.8;
        
        if any(isSrc)
            [x,y] = nxys(N,isSrc); [x,y]=applyOffset(x,y,off);
            scatter(x,y,gSize(isSrc),'s','filled','MarkerFaceColor',SourceFaceColor,'MarkerFaceAlpha',0.95);
        end
        if any(isDem)
            [x,y] = nxys(N,isDem); [x,y]=applyOffset(x,y,off);
            scatter(x,y,dSize(isDem),'o','filled','MarkerFaceColor',DemandFaceColor,'MarkerFaceAlpha',0.95);
        end
        if any(isOth)
            [x,y] = nxys(N,isOth); [x,y]=applyOffset(x,y,off);
            scatter(x,y,oSize(isOth),'^','filled','MarkerFaceColor',NetworkColor,'MarkerFaceAlpha',0.95);
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
        if isempty(CLim)
            mLims = pickLims([EM(:); NM(:)]);
        else
            mLims = CLim;
        end
        % --- EDGES (color by EdgeMetric) ---
        for e = 1:numel(S.Edge)
            [ex,ey] = getXY(S.Edge(e)); if isempty(ex), continue; end
            col = value2color(EM(e), mLims, flowmap, NetworkColor);
            [ex,ey] = applyOffset(ex,ey,off);
            plot(ex,ey,'-','Color',col,'LineWidth',EdgeWidth);
        end
        % --- NODES (color & size by NodeMetric) ---
        idx = true(1,numel(S.Node));
        [x_all,y_all] = nxys(S.Node, idx);
        [x,y] = applyOffset(x_all,y_all,off);
        C = zeros(numel(S.Node),3);
        for i=1:numel(S.Node)
            C(i,:) = value2color(NM(i), mLims, flowmap, NetworkColor);
        end
        sSize = mapToSize(NM, mLims, NodeMinSize, NodeMaxSize);
        scatter(ax, x, y, sSize, C, 'o', 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.95);
    end

    function plotLinks(field, fromK, offFrom, toK, offTo)
        if ~isfield(System, field) || isempty(System.(field)), return; end
        L = System.(field);
        if ~strcmp(Mode,'Topology')
            switch field
                case 'GasToPower'
                    if (strcmp(Mode,'Flow')&&strcmp(FlowSource,'custom'))
                        ipsVals = Interdependency.(field)(:,6)./Interdependency.(field)(:,5);
                    elseif strcmp(Mode,'Metric')
                        ipsVals = Interdependency.(field)(:,6);
                    else
                        ipsVals = L(:,5)./L(:,5);
                    end
                    
                otherwise
                    if (strcmp(Mode,'Flow')&&strcmp(FlowSource,'custom'))
                        ipsVals = Interdependency.(field)(:,4)./Interdependency.(field)(:,3);
                    elseif strcmp(Mode,'Metric')
                        ipsVals = Interdependency.(field)(:,4);
                    else
                        ipsVals = L(:,3)./L(:,3);
                    end
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
            plot([fx tx],[fy ty],'-','Color',col,'LineWidth',EdgeWidth);
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
        if idx>0 && idx<=numel(N) && hasXY(N(idx))
            x = N(idx).Longitude; y = N(idx).Latitude; [x,y] = applyOffset(x,y,off);
        else
            x = NaN; y = NaN;
        end
    end
end


function [xmin,xmax,ymin,ymax] = localBounds(Boundary, System)
% Compute overall bounding box of boundary and nodes
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
% resolveColormap Validate a user colormap (Mx3), clamp to [0,1], and resample to 256.
% If C is empty, return the fallback map. Always returns 256¡Á3.
if isempty(C)
    M = fallback; return;
end
if ~(isnumeric(C) && size(C,2)==3)
    error('Invalid %s: must be an M¡Á3 numeric matrix with values in [0,1].', name);
end
if any(~isfinite(C(:)))
    error('%s contains non-finite values.', name);
end
C = max(0,min(1,C));
if size(C,1) < 256
    M = resampleColormap(C,256);
else
    M = C;
end
end

function cmap2 = resampleColormap(cmap, targetN)
m = size(cmap,1);
xi = linspace(1,m,m);
xo = linspace(1,m,targetN);
cmap2 = [interp1(xi,cmap(:,1),xo,'linear').', ...
    interp1(xi,cmap(:,2),xo,'linear').', ...
    interp1(xi,cmap(:,3),xo,'linear').'];
cmap2 = max(0,min(1,cmap2));
end


function v = gopt(s, name, default)
if isfield(s,name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

function [ex,ey] = getXY(edge)
ex = []; ey = [];
if isfield(edge,'X') && isfield(edge,'Y')
    ex = edge.X; ey = edge.Y;
    if ~isnumeric(ex) || ~isnumeric(ey) || numel(ex)~=numel(ey)
        ex=[]; ey=[]; return;
    end
end
end

function tf = hasXY(node)
tf = isfield(node,'Longitude') && isfield(node,'Latitude') && ...
    isfinite(node.Longitude) && isfinite(node.Latitude);
end

function val = fget(s, field, default)
if isfield(s, field) && ~isempty(s.(field)) && isnumeric(s.(field))
    val = full(s.(field)); else, val = default; end
end

function c = value2color(v, lims, cmap, fallback)
if ~isfinite(v) || any(~isfinite(lims)) || lims(2)<=lims(1)
    c = fallback; return;
end
v = max(min(v, lims(2)), lims(1));
t = (v - lims(1)) / (lims(2) - lims(1));
k = 1 + floor(t*(size(cmap,1)-1));
k = max(1, min(k, size(cmap,1)));
c = cmap(k,:);
end

function sz = mapToSize(vals, lims, smin, smax)
vals(~isfinite(vals)) = NaN;
if isempty(lims) || any(~isfinite(lims)) || lims(2)<=lims(1)
    lims = [0 1];
end
vals = abs(vals); % size uses absolute value
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
    w = abs(w);
    lo = min(w); hi = max(w);
    if hi <= lo, hi = lo + eps; end
    lims = [lo hi];
end
end

function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end

function [x,y] = nxys(nodes,mask)
idx = find(mask);
x = arrayfun(@(i) nodes(i).Longitude, idx);
y = arrayfun(@(i) nodes(i).Latitude,  idx);
end

function v = getOptVec(s, name, n, defaultVal, label)
if isfield(s,name) && ~isempty(s.(name))
    v0 = s.(name)(:);
    if numel(v0) ~= n
        error('%s length (%d) does not match expected length (%d).', label, numel(v0), n);
    else
        v = v0;
    end
else
    warning('options.%s is missing. Filling with %d.',label, defaultVal);
    v = repmat(defaultVal, n, 1);
end
end
