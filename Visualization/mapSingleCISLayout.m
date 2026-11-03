function mapSingleCISLayout(Boundary, CIS, Mode, options)
% INTRODUCTION:
%   mapSingleCISLayout renders one CIS within an optional geographic
%   boundary. It supports three visualization modes:
%     - 'Topology' - structural view of nodes and edges.
%     - 'Flow'     - edge flow colors and node generation/demand sizes.
%     - 'Metric'   - user-provided edge and node metrics.
%
% INPUTS:
%   Boundary : scalar struct with optional X and Y coordinate vectors. If
%              both fields are present and nonempty, the boundary is drawn.
%   CIS      : struct with CIS.Node(1-by-N) and CIS.Edge(1-by-E). Nodes
%              require Longitude and Latitude; edges require X and Y.
%   Mode     : 'Topology' | 'Flow' | 'Metric'
%   options  : struct controlling styles, data sources, and metrics:
%     General (all modes)
%       - SystemType    : 'power' | 'gas' | 'water' | 'road' (required)
%       - FigureColor   : figure background color (default 'white'); used
%                         only when AxesOpt is not supplied
%       - BoundaryColor : boundary RGB color (default [0 0 0])
%       - NetworkColor  : neutral network RGB color (default [96 96 96]/255)
%       - EdgeWidth     : network edge width (default 1)
%       - AxesOpt       : handle to existing axes for drawing
%
%     Topology mode
%       - MarkerSizeTopo : node marker size (default 4)
%
%     Topology and Metric modes
%       - Topology : 'nodes' | 'edges' | 'both' (default 'both')
%
%     Flow and Metric modes
%       - FlowColormap : nonempty M-by-3 colormap in [0,1]; despite its
%                        historical name, it is also used in Metric mode
%       - NodeMinSize/NodeMaxSize : scatter size range (default 10/20)
%       - ValueLabel   : colorbar label (default '|Flow|' or 'Metric')
%
%     Flow mode
%       - FlowSource : 'system' (default) | 'custom'; road systems do not
%                      support Flow mode
%       - EdgeFlow   : E-by-1 custom edge-flow vector
%       - NodeGen    : N-by-1 custom generation vector
%       - NodeDem    : N-by-1 custom demand vector
%       - SourceFaceColor : source-node RGB color (default [0.15 0.47 0.80])
%       - DemandFaceColor : demand-node RGB color (default [0.86 0.33 0.10])
%         System Flow uses RealFlow, RealGeneration, and RealDemand. A node
%         with both generation and demand is displayed as a source node.
%
%     Metric mode
%       - EdgeMetric   : E-by-1 numeric vector (optional)
%       - NodeMetric   : N-by-1 numeric vector (optional)
%       - MetricLimits : [min max] color limits (default finite metric range)
%         Missing or nonfinite values use NetworkColor and minimum size;
%         they are not interpreted as zero-valued metrics.
%
% OUTPUT:
%   The function draws on options.AxesOpt when supplied; otherwise it creates
%   a new figure and axes. No variables are returned.

%% Step 1: required args & basic checks
if nargin < 3 || isempty(Mode)
    error('Argument "mode" is required and must be "Topology", "Flow" or "Metric".');
end
if nargin < 4 || isempty(options), options = struct(); end

Mode = string(Mode);
validModes = ["Topology","Flow","Metric"];
if ~any(Mode == validModes)
    error('Invalid mode: %s. Use "Topology", "Flow" or "Metric".', Mode);
end

%% Step 2: Validate inputs
SystemType = string(gopt(options,'SystemType',''));
if strlength(SystemType) == 0
    error('options.SystemType is required: power | gas | water | road.');
end
validTypes = ["power","gas","water","road"];
if ~any(SystemType == validTypes)
    error('options.SystemType must be one of: power, gas, water, road.');
end

if SystemType == "road" && Mode == "Flow"
    error('SystemType="road" only supports Topology and Metric mode.');
end

FigureColor     = gopt(options,'FigureColor','white');
BoundaryColor   = validateRGB(gopt(options,'BoundaryColor',[0 0 0]),'BoundaryColor');
NetworkColor    = validateRGB(gopt(options,'NetworkColor',[96 96 96]/255),'NetworkColor');
MarkerSizeTopo  = gopt(options,'MarkerSizeTopo',4);
EdgeWidth       = gopt(options,'EdgeWidth',1);
AxesOpt         = gopt(options,'AxesOpt',[]);
Topology        = lower(string(gopt(options,'Topology','both')));
FlowColormap    = gopt(options,'FlowColormap',[]);
NodeMinSize     = gopt(options,'NodeMinSize',10);
NodeMaxSize     = gopt(options,'NodeMaxSize',20);
ValueLabel      = gopt(options,'ValueLabel','');
SourceFaceColor = gopt(options,'SourceFaceColor',[0.15 0.47 0.80]);
DemandFaceColor = gopt(options,'DemandFaceColor',[0.86 0.33 0.10]);
CLim            = gopt(options,'MetricLimits',[]);

if ~isnumeric(EdgeWidth) || ~isscalar(EdgeWidth) || ~isfinite(EdgeWidth) || EdgeWidth <= 0
    error('options.EdgeWidth must be a positive finite scalar.');
end
if any(Mode == ["Topology","Metric"]) && (~isscalar(Topology) || ~any(Topology == ["nodes","edges","both"]))
    error('options.Topology must be "nodes", "edges" or "both".');
end
if Mode == "Topology" && (~isnumeric(MarkerSizeTopo) || ~isscalar(MarkerSizeTopo) || ...
        ~isfinite(MarkerSizeTopo) || MarkerSizeTopo <= 0)
    error('options.MarkerSizeTopo must be a positive finite scalar.');
end
if any(Mode == ["Flow","Metric"])
    if ~(ischar(ValueLabel) || (isstring(ValueLabel) && isscalar(ValueLabel)))
        error('options.ValueLabel must be a character vector or string scalar.');
    end
    if ~isnumeric(NodeMinSize) || ~isscalar(NodeMinSize) || ~isfinite(NodeMinSize) || NodeMinSize <= 0 || ...
            ~isnumeric(NodeMaxSize) || ~isscalar(NodeMaxSize) || ~isfinite(NodeMaxSize) || NodeMaxSize < NodeMinSize
        error('options.NodeMinSize and NodeMaxSize must be positive finite scalars with NodeMaxSize >= NodeMinSize.');
    end
end
if Mode == "Flow"
    SourceFaceColor = validateRGB(SourceFaceColor,'SourceFaceColor');
    DemandFaceColor = validateRGB(DemandFaceColor,'DemandFaceColor');
end

if Mode ~= "Topology"
    if ~isempty(FlowColormap)
        if ~isnumeric(FlowColormap) || ~ismatrix(FlowColormap) || size(FlowColormap,2) ~= 3 || ...
                ~isreal(FlowColormap) || any(~isfinite(FlowColormap(:))) || ...
                any(FlowColormap(:) < 0 | FlowColormap(:) > 1)
            error('options.FlowColormap must be a nonempty M-by-3 numeric matrix with values in [0,1].');
        end
        if size(FlowColormap,1) < 256
            cmap = resampleColormap(FlowColormap,256);
        else
            cmap = FlowColormap;
        end
    elseif Mode == "Flow"
        cmap = parula(256);
    else
        cmap = autumn(256);
    end
end

FlowSource = string(gopt(options,'FlowSource','system'));
if Mode == "Flow"
    if ~isscalar(FlowSource) || ~any(FlowSource == ["system","custom"])
        error('options.FlowSource must be "system" or "custom" when mode="Flow".');
    end
end

N = numel(CIS.Node); E = numel(CIS.Edge);
if Mode == "Flow" && FlowSource == "custom"
    EdgeFlow = getOptVec(options,'EdgeFlow',E,'EdgeFlow');
    NodeGen  = getOptVec(options,'NodeGen',N,'NodeGen');
    NodeDem  = getOptVec(options,'NodeDem',N,'NodeDem');
else
    EdgeFlow = []; NodeGen = []; NodeDem = [];
end

if Mode == "Metric"
    if isfield(options,'EdgeMetric') && ~isempty(options.EdgeMetric)
        EdgeMetric = getOptVec(options,'EdgeMetric',E,'EdgeMetric');
    else
        EdgeMetric = NaN(E,1);
    end
    if isfield(options,'NodeMetric') && ~isempty(options.NodeMetric)
        NodeMetric = getOptVec(options,'NodeMetric',N,'NodeMetric');
    else
        NodeMetric = NaN(N,1);
    end
    if ~isempty(CLim)
        if ~isnumeric(CLim) || numel(CLim) ~= 2 || ~isreal(CLim) || any(~isfinite(CLim(:)))
            error('options.MetricLimits must be a finite 1-by-2 numeric vector [min max].');
        end
        CLim = CLim(:).';
        if CLim(1) >= CLim(2)
            error('options.MetricLimits must satisfy MetricLimits(1) < MetricLimits(2).');
        end
    end
end

%% Step 3: Figure & Axes
if isempty(AxesOpt)
    figure('Color',FigureColor);
    ax = axes;
else
    ax = AxesOpt;
end
holdState = ishold(ax);
hold(ax,'on'); grid(ax,'on');

if isstruct(Boundary) && isscalar(Boundary) && isfield(Boundary,'X') && isfield(Boundary,'Y') && ...
        ~isempty(Boundary.X) && ~isempty(Boundary.Y)
    plot(ax, Boundary.X, Boundary.Y, '-', 'Color', BoundaryColor, 'LineWidth', 1.5);
end

%% Step 4: Plot by mode
switch Mode
    case "Topology"
        if any(Topology == ["edges","both"])
            for e = 1:E
                [ex,ey] = getXY(CIS.Edge(e)); if isempty(ex), continue; end
                plot(ax, ex, ey, '-', 'Color', NetworkColor, 'LineWidth', EdgeWidth);
            end
        end
        if any(Topology == ["nodes","both"])
            for n = 1:N
                td = fget(CIS.Node(n),'TargetDemand',0);
                markerStyle = iff(td > 0, 'o', 's');
                plot(ax, CIS.Node(n).Longitude, CIS.Node(n).Latitude, markerStyle, ...
                    'MarkerSize', MarkerSizeTopo, ...
                    'Color', NetworkColor, 'MarkerFaceColor', NetworkColor);
            end
        end
    case "Flow"
        if FlowSource == "custom"
            flowVals = EdgeFlow(:);
            genVals  = NodeGen(:);
            demVals  = NodeDem(:);
        else
            flowVals = arrayfun(@(e) fget(e,'RealFlow',NaN), CIS.Edge);
            genVals  = arrayfun(@(n) fget(n,'RealGeneration',0), CIS.Node);
            demVals  = arrayfun(@(n) fget(n,'RealDemand',0), CIS.Node);
        end
        
        flowForColor = abs(flowVals);

        flowLims = pickLims(flowForColor);
        genLims  = pickLims(genVals);
        demLims  = pickLims(demVals);

        for e = 1:E
            [ex,ey] = getXY(CIS.Edge(e)); if isempty(ex), continue; end
            c  = value2color(flowForColor(e), flowLims, cmap, NetworkColor);
            plot(ax, ex, ey, '-', 'Color', c, 'LineWidth', EdgeWidth);
        end
        colormap(ax, cmap);
        cb = colorbar('peer', ax);
        if isempty(ValueLabel), ValueLabel = '|Flow|'; end
        cb.Label.String = ValueLabel;
        cb.Label.FontSize = 14;
        cb.FontSize = 14;
        caxis(ax, flowLims);

        % Source role takes precedence when a node also has demand.
        isSource = genVals > 0;
        isDemand = ~isSource & (demVals > 0);
        isOther  = ~isSource & ~isDemand;
        
        sSize = mapToSize(genVals, genLims, NodeMinSize, NodeMaxSize);
        dSize = mapToSize(demVals, demLims, NodeMinSize, NodeMaxSize);

        if any(isSource)
            [x,y] = nxys(CIS.Node,isSource);
            scatter(ax, x, y, sSize(isSource), 's', 'filled', ...
                'MarkerFaceColor', SourceFaceColor, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','none');
        end
        if any(isDemand)
            [x,y] = nxys(CIS.Node,isDemand);
            scatter(ax, x, y, dSize(isDemand), 'o', 'filled', ...
                'MarkerFaceColor', DemandFaceColor, 'MarkerFaceAlpha', 0.95, 'MarkerEdgeColor','none');
        end
        if any(isOther)
            [x,y] = nxys(CIS.Node,isOther);
            scatter(ax, x, y, NodeMinSize*0.8, '^', 'filled', ...
                'MarkerFaceColor', NetworkColor, 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor','none');
        end
        
    case "Metric"
        allVals = [EdgeMetric(:); NodeMetric(:)];
        if isempty(CLim)
            mLims   = pickLims(allVals);
        else
            mLims = CLim;
        end
        
        if any(Topology == ["edges","both"])
            for e = 1:E
                [ex,ey] = getXY(CIS.Edge(e)); if isempty(ex), continue; end
                col = value2color(EdgeMetric(e), mLims, cmap, NetworkColor);
                plot(ax, ex, ey, '-', 'Color', col, 'LineWidth', EdgeWidth);
            end
        end
        
        if any(Topology == ["nodes","both"])
            [x_all,y_all] = nxys(CIS.Node, true(1,N));
            C = zeros(N,3);
            for i=1:N
                C(i,:) = value2color(NodeMetric(i), mLims, cmap, NetworkColor);
            end
            sSize = mapToSize(NodeMetric, mLims, NodeMinSize, NodeMaxSize);
            scatter(ax, x_all, y_all, sSize, C, 'o', 'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.95);
        end
        colormap(ax, cmap);
        cb = colorbar('peer', ax);
        if isempty(ValueLabel), ValueLabel = 'Metric'; end
        cb.Label.String = ValueLabel;
        cb.Label.FontSize = 16;
        cb.FontSize = 16;
        cb.FontName = 'Times New Roman';
        caxis(ax, mLims);
end

%% Step 5: View
axis(ax,'equal');
axis(ax,'off'); 
axis(ax,'tight'); 
if ~holdState, hold(ax,'off'); end

end



function cmap2 = resampleColormap(cmap, targetN)
m = size(cmap,1);
if m == 1
    cmap2 = repmat(cmap,targetN,1);
    return;
end
xo = linspace(1,m,targetN);
cmap2 = [interp1(1:m,cmap(:,1),xo,'linear').', ...
    interp1(1:m,cmap(:,2),xo,'linear').', ...
    interp1(1:m,cmap(:,3),xo,'linear').'];
end

function v = gopt(s, name, default)
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end

function c = validateRGB(c, name)
if ~isnumeric(c) || numel(c) ~= 3 || ~isreal(c) || any(~isfinite(c(:))) || any(c(:) < 0 | c(:) > 1)
    error('options.%s must be a finite RGB triplet with values in [0,1].',name);
end
c = reshape(c,1,3);
end

function [ex,ey] = getXY(edge)
ex = edge.X;
ey = edge.Y;
end

function val = fget(s, field, default)
if isfield(s, field) && ~isempty(s.(field)) && isnumeric(s.(field))
    val = full(s.(field));
else
    val = default;
end
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
invalid = ~isfinite(vals);
sz = (vals - lims(1)) / (lims(2)-lims(1));
sz = max(0, min(1, sz));
sz = smin + sz .* (smax - smin);
sz(invalid) = smin;
end

function lims = pickLims(v)
w = v(isfinite(v));
if isempty(w)
    lims = [0 1];
else
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
x = [nodes(mask).Longitude];
y = [nodes(mask).Latitude];
end

function v = getOptVec(s, name, n, label)
if isfield(s,name) && ~isempty(s.(name))
    v0 = s.(name)(:);
    if ~isnumeric(v0) || ~isreal(v0)
        error('options.%s must be a real numeric vector.',label);
    elseif numel(v0) ~= n
        error('%s length (%d) does not match expected length (%d).',label,numel(v0),n);
    else
        v = v0;
    end
else
    error('options.%s is required and must not be empty.',label);
end
end

