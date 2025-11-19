function mapTerminalZoneValue(TerminalZone, Value, options)
% INTRODUCTION
%   Draws colored polygons from a TerminalZone struct using fill(). Each
%   cell is shaded according to its corresponding Value entry. Designed for
%   1¡ÁN TerminalZone arrays (fields .X, .Y). NaN values use NaNColor. If
%   options.AxesOpt is given, drawing occurs on that axes; otherwise a new
%   figure is created. A colorbar labeled by ValueLabel is added
%   automatically.

% INPUT
%   TerminalZone : 1¡ÁN struct array. Each element must include:
%                  .X (1¡ÁK double) ¡ª polygon vertex X coordinates (may end with NaN)
%                  .Y (1¡ÁK double) ¡ª polygon vertex Y coordinates (may end with NaN)
%                  Additional fields (if any) are ignored.
%
%   Value        : N¡Á1 or 1¡ÁN numeric vector, one scalar per polygon.
%                  Finite entries are used for color mapping; NaN entries use NaNColor.
%
%   options      : (optional) struct controlling appearance and behavior:
%       - FigureColor : figure background color (default = 'white')
%       - ValueLimits : [vmin vmax] (default = finite data range)
%       - Colormap    : M¡Á3 colormap matrix (default = parula(256))
%       - NaNColor    : 1¡Á3 RGB triplet for NaN polygons (default = [0.8 0.8 0.8])
%       - EdgeColor   : RGB triplet for polygon edges (default = [0.4 0.4 0.4])
%       - EdgeAlpha   : edge transparency, scalar in [0,1] (default = 0.4)
%       - FaceAlpha   : face transparency, scalar in [0,1] (default = 1)
%       - AxesOpt     : axes handle; if provided, draw on this axes instead of creating a new figure
%       - ValueLabel  : string used as colorbar label (default = 'Value')
% 
% OUTPUT
%   none ¡ª this is a renderer that draws directly to a figure or provided axes.
%          If AxesOpt is not provided, a new figure is created. A colorbar
%          is added with label specified by ValueLabel, and axes are set to
%          equal and tight with axis visibility turned off.

%% Step 1: basic checks
if nargin < 3, options = struct(); end
N = numel(TerminalZone);
if numel(Value) ~= N
    error('Length of Value (%d) must match number of TerminalZone cells (%d).', numel(Value), N);
end
Value = Value(:);

%% Step 2: parse options
FigureColor  = gopt(options,'FigureColor','white');
CLim         = gopt(options,'ValueLimits', []);
Cmap         = gopt(options,'Colormap', parula(256));
NaNColor     = gopt(options,'NaNColor', [0.8 0.8 0.8]);
EdgeColor    = gopt(options,'EdgeColor',[0.4 0.4 0.4]);
EdgeAlpha    = gopt(options,'EdgeAlpha',0.4);
FaceAlpha    = gopt(options,'FaceAlpha',1);
AxesOpt      = gopt(options,'AxesOpt',[]);
ValueLabel   = gopt(options,'ValueLabel','Value'); 

%% Step 3: figure/axes setup
if isempty(AxesOpt)
    figure('Color',FigureColor);
    ax = axes; hold(ax,'on');
else
    ax = AxesOpt; axes(ax); hold(ax,'on');
end

if isempty(CLim)
    finiteVals = Value(isfinite(Value));
    if isempty(finiteVals)
        CLim = [0 1];
    else
        CLim = [min(finiteVals) max(finiteVals)];
        if CLim(1)==CLim(2), CLim = CLim + [-0.5 0.5]; end
    end
else
    if ~isnumeric(CLim) || numel(CLim) ~= 2
        error('options.ValueLimits must be a 1¡Á2 numeric vector [min max].');
    end
    if CLim(1) >= CLim(2)
        error('options.ValueLimits must satisfy ValueLimits(1) < ValueLimits(2).');
    end
end
if size(Cmap,1) < 256
    Cmap = resampleColormap(Cmap,256);
end
colormap(ax,Cmap);
caxis(ax,CLim);
%% Step 4: draw polygons
nCmap = size(Cmap,1);
for k = 1:N
    if isempty(TerminalZone(k).X) || isempty(TerminalZone(k).Y), continue; end
    v = Value(k);
    if isnan(v)
        rgb = NaNColor;
    else
        idx = 1 + round((nCmap-1)*(v-CLim(1))/max(eps,CLim(2)-CLim(1)));
        idx = max(1,min(nCmap,idx));
        rgb = Cmap(idx,:);
    end

    X = TerminalZone(k).X; Y = TerminalZone(k).Y;
    if numel(X)>1 && isnan(X(end)) && isnan(Y(end))
        X = X(1:end-1); Y = Y(1:end-1);
    end
        fill(ax, X, Y,  rgb, 'EdgeColor',EdgeColor, 'EdgeAlpha',EdgeAlpha, ...
             'FaceAlpha',FaceAlpha, 'LineWidth',0.8);
end

axis(ax,'equal');
axis(ax,'tight'); box(ax,'on');
colorbar(ax); 
ylabel(colorbar,ValueLabel);
hold(ax,'off');
axis(ax,'off');
end

function out = gopt(s,name,default)
if isfield(s,name) && ~isempty(s.(name)), out = s.(name); else, out = default; end
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
