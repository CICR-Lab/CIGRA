function mapSingleCISLayout(Boundary, CIS, Mode, options)
% INTRODUCTION:
%   mapSingleCISLayout visualizes a single Critical Infrastructure System (CIS) ¡ª such as
%   power, gas, water, or road networks ¡ª within a geographic boundary. The function
%   provides three visualization modes to support different analytical purposes:
%
%     1. 'Topology' ¡ª Displays only the network structure, i.e., nodes and edges.
%     2. 'Flow'     ¡ª Depicts flow magnitudes (e.g., power, gas, or water) along edges
%                     and node generation/demand. Flow values can come from system data
%                     or user-defined custom inputs.
%     3. 'Metric'   ¡ª Visualizes user-defined scalar metrics (e.g., component failure
%                     probability, importance, or 0¨C1 operational states) for nodes and
%                     edges using color gradients.
%
% INPUTS:
%   Boundary : struct with optional fields X, Y ¡ª coordinates outlining the study area.
%              If provided, the boundary is drawn as a thin outline in the background.
%
%   CIS      : struct describing the system topology with two arrays:
%              - CIS.Node(i): must include Longitude, Latitude, and optionally fields such as
%                  ID, TargetDemand, MaxGeneration, RealDemand, RealGeneration.
%              - CIS.Edge(e): must include FromNodeID, ToNodeID, X, Y (polyline coordinates),
%                  and optionally RealFlow, Capacity, etc.
%
%   Mode     : string specifying visualization type:
%              - 'Topology' ¡ª structural view of network.
%              - 'Flow'     ¡ª edge flows and node generation/demand.
%              - 'Metric'   ¡ª user-defined per-node and per-edge scalar metrics.
%
%   options  : struct of configuration parameters. The following fields are supported:
%
%       Required:
%         - SystemType      : 'power' | 'gas' | 'water' | 'road'
%
%       General (all modes):
%         - FigureColor     : background color (default 'white')
%         - BoundaryColor   : boundary line color (default [0 0 0])
%         - NetworkColor    : base color for network edges (default [96 96 96]/255)
%         - MarkerSizeTopo  : node marker size for topology view (default 4)
%         - EdgeWidthTopo   : line width for topology view (default 1)
%         - AxesOpt         : axes handle; if provided, draw on this axes instead of creating a new figure
%         - FlowColormap    : M¡Á3 colormap matrix in [0,1], optional.
%
%       Flow mode only:
%         - FlowSource      : 'system' (use CIS data) | 'custom' (use user inputs)
%         - EdgeFlow        : E¡Á1 numeric vector of custom edge flow values (required if
%                             FlowSource='custom')
%         - NodeGen, NodeDem: N¡Á1 numeric vectors of custom generation/demand values.
%         - SourceFaceColor : RGB color for generation/source nodes (default [0.15 0.47 0.80])
%         - DemandFaceColor : RGB color for demand/sink nodes (default [0.86 0.33 0.10])
%         - NodeMinSize / NodeMaxSize : scatter size scaling range.
%
%       Metric mode only:
%         - EdgeMetric, NodeMetric : REQUIRED numeric vectors (E¡Á1 and N¡Á1, respectively).
%                                   Values can represent probabilities, performance indices,
%                                   or binary states (0/1). All are mapped to colors.
%         - EdgeWidthMetric  : line width for edges in metric mode (default 1.2)
%         - NodeMinSize / NodeMaxSize : scatter size scaling range.
% 
%       Topology mode and Metric mode:
%         - Topology        : select from 'nodes','edges','both'.
%
% OUTPUT:
%   The function produces a figure displaying the CIS network in the specified mode.
%   It returns no explicit variables but renders plots with:
%     - Colored edges and nodes according to flow magnitudes or metric values.
%     - Automatically scaled colorbars and axis limits.
%     - Distinct legends or markers for different node types in Flow mode.

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
% Validate SystemType
SystemType = string(gopt(options,'SystemType',''));
if strlength(SystemType) == 0
    error('options.SystemType is required: power | gas | water | road.');
end
validTypes = ["power","gas","water","road"];
if ~any(SystemType == validTypes)
    error('options.SystemType must be one of: power, gas, water, road.');
end

% Enforce road Mode Topology
if SystemType == "road" && Mode == "Flow"
    error('SystemType="road" only supports Topology and Metric mode.');
end

% Load style / behavior options
FigureColor     = gopt(options,'FigureColor','white');
BoundaryColor   = gopt(options,'BoundaryColor',[0 0 0]);
NetworkColor    = gopt(options,'NetworkColor',[96 96 96]/255);
MarkerSizeTopo  = gopt(options,'MarkerSizeTopo',4);
EdgeWidth       = gopt(options,'EdgeWidth',1);
AxesOpt         = gopt(options,'AxesOpt',[]);
Topology        = gopt(options,'Topology', 'both'); % nodes and edges

% When choose Mode=Flow
FlowColormap    = gopt(options,'FlowColormap',[]);
NodeMinSize     = gopt(options,'NodeMinSize',10);
NodeMaxSize     = gopt(options,'NodeMaxSize',20);
SourceFaceColor = gopt(options,'SourceFaceColor',[0.15 0.47 0.80]);
DemandFaceColor = gopt(options,'DemandFaceColor',[0.86 0.33 0.10]);

% When choose Mode=Metric
CLim         = gopt(options,'MetricLimits', []);

% Colormap for edge flow
if ~isempty(FlowColormap)
    if ~(isnumeric(FlowColormap) && size(FlowColormap,2)==3)
        error('Invalid FlowColormap: must be an M¡Á3 numeric matrix with values in [0,1].');
    end
    if any(~isfinite(FlowColormap(:)))
        error('FlowColormap contains non-finite values.');
    end
    if min(FlowColormap(:))<0 || max(FlowColormap(:))>1
        warning('FlowColormap values out of [0,1]. Clamping automatically.');
        FlowColormap = max(0,min(1,FlowColormap));
    end
    if size(FlowColormap,1) < 256
        cmap = resampleColormap(FlowColormap,256);
    else
        cmap = FlowColormap;
    end
else
    switch Mode
        case 'Flow'
            cmap = parula(256);
        case 'Metric'
            cmap = autumn(256);
    end
end

% FlowSource handling: only required/validated when mode=="Flow"
FlowSourceRaw = gopt(options,'FlowSource','');
if strcmp(Mode,'Flow')
    if isempty(FlowSourceRaw)
        warning('mode="Flow" without options.FlowSource. Defaulting to "system" (system-provided values).');
        FlowSource = "system";
    else
        FlowSource = string(FlowSourceRaw);
        if ~any(FlowSource == ["system","custom"])
            error('options.FlowSource must be "system" or "custom" when mode="Flow".');
        end
    end
else
    FlowSource = NaN;
end

% Validate custom inputs if requested
N = numel(CIS.Node); E = numel(CIS.Edge);
if strcmp(FlowSource,'custom') && strcmp(Mode,'Flow')
    % Custom vectors; pad/trunc as needed with warnings
    EdgeFlow = getOptVec(options,'EdgeFlow',E, 0, 'EdgeFlow');
    NodeGen  = getOptVec(options,'NodeGen', N, 0, 'NodeGen');
    NodeDem  = getOptVec(options,'NodeDem', N, 0, 'NodeDem');
else
    EdgeFlow = []; NodeGen = []; NodeDem = [];
end

% When choose Mode=Metric
if strcmp(Mode,'Metric')
    EdgeMetric = getOptVec(options,'EdgeMetric',E, 0, 'EdgeMetric');
    NodeMetric = getOptVec(options,'NodeMetric',N, 0, 'NodeMetric');
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
if isempty(AxesOpt)
    figure('Color',FigureColor);
    ax = axes; hold(ax,'on');
else
    ax = AxesOpt; axes(ax); hold(ax,'on');
end
hold(ax,'on'); grid(ax,'on');

% Boundary
if ~isempty(Boundary.X)
    plot(ax, Boundary.X, Boundary.Y, '-', 'Color', BoundaryColor, 'LineWidth', 1);
end

%% Step 4: Plot by mode
switch Mode
    case "Topology"
        % --- Edges ---
        if ismember(Topology,{'edges','both'})
            for e = 1:E
                [ex,ey] = getXY(CIS.Edge(e)); if isempty(ex), continue; end
                plot(ax, ex, ey, '-', 'Color', NetworkColor, 'LineWidth', EdgeWidth);
            end
        end
        % --- Nodes ---
        if ismember(Topology,{'nodes','both'})
            for n = 1:N
                td = fget(CIS.Node(n),'TargetDemand',0);
                markerStyle = iff(td > 0, 'o', 's');
                plot(ax, CIS.Node(n).Longitude, CIS.Node(n).Latitude, markerStyle, ...
                    'MarkerSize', MarkerSizeTopo, ...
                    'Color', NetworkColor, 'MarkerFaceColor', NetworkColor);
            end
        end
    case "Flow"
        % --- Gather values (system or custom) ---
        if FlowSource == "custom"
            flowVals = EdgeFlow(:);
            genVals  = NodeGen(:);
            demVals  = NodeDem(:);
        else
            flowVals = arrayfun(@(e) fget(e,'RealFlow',NaN), CIS.Edge);
            genVals  = arrayfun(@(n) fget(n,'MaxGeneration',0), CIS.Node);
            demVals  = arrayfun(@(n) fget(n,'TargetDemand',0), CIS.Node);
        end
        
        % --- Normalize for color ---
        flowForColor = abs(flowVals);
        flowForColor(~isfinite(flowForColor)) = NaN;
        
        % --- Limits  ---
        flowLims = pickLims(flowForColor);
        genLims  = pickLims(genVals);
        demLims  = pickLims(demVals);
        
        % --- EDGES (color by flowForColor) ---
        for e = 1:E
            [ex,ey] = getXY(CIS.Edge(e)); if isempty(ex), continue; end
            c  = value2color(flowForColor(e), flowLims, cmap, NetworkColor);
            plot(ax, ex, ey, '-', 'Color', c, 'LineWidth', EdgeWidth);
        end
        colormap(ax, cmap);
        cb = colorbar('peer', ax);
        cb.Label.String = '|Flow|';
        if all(isfinite(flowLims)) && diff(flowLims) > 0
            caxis(ax, flowLims);
        end
        
        % --- NODES (size & color by type/value) ---
        isSource = genVals > 0;
        isDemand = ~isSource & (demVals > 0);
        isOther  = ~isSource & ~isDemand;
        
        sSize = mapToSize(genVals, genLims, NodeMinSize, NodeMaxSize);
        dSize = mapToSize(demVals, demLims, NodeMinSize, NodeMaxSize);
        oSize = ones(size(genVals))*NodeMinSize*0.8;
        
        % sources
        if any(isSource)
            [x,y] = nxys(CIS.Node,isSource);
            scatter(ax, x, y, sSize(isSource), 's', 'filled', ...
                'MarkerFaceColor', SourceFaceColor, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','none');
        end
        % demands
        if any(isDemand)
            [x,y] = nxys(CIS.Node,isDemand);
            scatter(ax, x, y, dSize(isDemand), 'o', 'filled', ...
                'MarkerFaceColor', DemandFaceColor, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','none');
        end
        % others
        if any(isOther)
            [x,y] = nxys(CIS.Node,isOther);
            scatter(ax, x, y, oSize(isOther), '^', 'filled', ...
                'MarkerFaceColor', NetworkColor, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor','none');
        end
        
    case "Metric"
        % Combined limits across node/edge to use a single colorbar
        allVals = [EdgeMetric(:); NodeMetric(:)];
        if isempty(CLim)
            mLims   = pickLims(allVals);
        else
            mLims = CLim;
        end
        
        % --- EDGES (color by EdgeMetric) ---
        if ismember(Topology,{'edges','both'})
            for e = 1:E
                [ex,ey] = getXY(CIS.Edge(e)); if isempty(ex), continue; end
                col = value2color(EdgeMetric(e), mLims, cmap, NetworkColor);
                plot(ax, ex, ey, '-', 'Color', col, 'LineWidth', EdgeWidth);
            end
        end
        
        % --- NODES (color & size by NodeMetric) ---
        if ismember(Topology,{'nodes','both'})
            idx = true(1,N);
            [x_all,y_all] = nxys(CIS.Node, idx);
            C = zeros(N,3);
            for i=1:N
                C(i,:) = value2color(NodeMetric(i), mLims, cmap, NetworkColor);
            end
            sSize = mapToSize(NodeMetric, mLims, NodeMinSize, NodeMaxSize);
            scatter(ax, x_all, y_all, sSize, C, 'o', 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.95);
        end
        colormap(ax, cmap);
        cb = colorbar('peer', ax);
        cb.Label.String = 'Metric';
        if all(isfinite(mLims)) && diff(mLims) > 0, caxis(ax, mLims); end
end

%% Step 5: View
axis(ax,'equal');
axis(ax,'off'); axis(ax,'tight'); hold(ax,'off');

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
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
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

function val = fget(s, field, default)
if isfield(s, field) && ~isempty(s.(field)) && isnumeric(s.(field))
    val = full(s.(field));
else
    val = default;
end
end

function c = value2color(v, lims, cmap, fallback)
% Map scalar v to RGB using colormap and limits; fallback if invalid
if ~isfinite(v) || isempty(lims) || any(~isfinite(lims)) || lims(2)<=lims(1)
    c = fallback; return;
end
% clamp
v = max(min(v, lims(2)), lims(1));
% normalize to [0,1]
t = (v - lims(1)) / (lims(2) - lims(1));
% index
k = 1 + floor(t*(size(cmap,1)-1));
k = max(1, min(k, size(cmap,1)));
c = cmap(k,:);
end

function sz = mapToSize(vals, lims, smin, smax)
vals(~isfinite(vals)) = NaN;
if isempty(lims) || any(~isfinite(lims)) || lims(2)<=lims(1)
    % fall back to percentiles
    v = vals; v = v(isfinite(v));
    if isempty(v)
        lims = [0 1];
    else
        lims = [prctile(v,5) prctile(v,95)];
        if lims(2)<=lims(1), lims = [min(v) max(v) + eps]; end
    end
end
% linear map into [smin, smax]
sz = (vals - lims(1)) / (lims(2)-lims(1));
sz = max(0, min(1, sz));
sz = smin + sz .* (smax - smin);
% replace NaNs
sz(~isfinite(sz)) = smin;
end

function lims = pickLims(v)
w = v(isfinite(v));
if isempty(w)
    lims = [0 1];
else
    if ~isvector(w), w = w(:); end
    lims = [min(w) max(w)];
    if lims(2) <= lims(1)
        lims(2) = lims(1) + eps;
    end
end
end

function out = iff(cond, a, b)
if cond, out = a; else, out = b; end
end

function [x,y] = nxys(nodes,mask)
idx = find(mask);
if isempty(idx)
    x = []; y = []; return;
end
x = arrayfun(@(i) nodes(i).Longitude, idx);
y = arrayfun(@(i) nodes(i).Latitude,  idx);
end

function v = getOptVec(s, name, n, defaultVal, label)
if isfield(s,name) && ~isempty(s.(name))
    v0 = s.(name)(:);
    if numel(v0) ~= n
        error('%s length (%d) does not match expected length (%d). Using best-effort trunc/pad.', label, numel(v0), n);
    else
        v = v0;
    end
else
    warning('options.%s is missing. Filling with %d.',label, defaultVal);
    v = repmat(defaultVal, n, 1);
end
end

