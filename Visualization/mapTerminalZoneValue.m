function mapTerminalZoneValue(TerminalZone, Value, options)
% INTRODUCTION
%   Draws colored polygons from a TerminalZone struct using fill(). Each
%   cell is shaded according to its corresponding Value entry. Designed for
%   1-by-N TerminalZone arrays (fields .X, .Y). If options.AxesOpt is
%   given, drawing occurs on that axes; otherwise a new figure is created.

% INPUT
%   TerminalZone : 1-by-N struct array. Each element must include:
%                  .X (1-by-K double) - polygon vertex X coordinates (may end with NaN)
%                  .Y (1-by-K double) - polygon vertex Y coordinates (may end with NaN)
%                  Additional fields (if any) are ignored.
%
%   Value        : N-by-1 or 1-by-N numeric vector, one scalar per polygon.
%                  Finite entries are used for color mapping; NaN and Inf use NaNColor.
%
%   options      : (optional) struct controlling appearance and behavior:
%       - FigureColor : figure background color (default = 'white')
%       - ValueLimits : [vmin vmax] (default = finite data range)
%       - Colormap    : M-by-3 colormap matrix (default = parula(256))
%       - NaNColor    : 1-by-3 RGB triplet for NaN polygons (default = [0.8 0.8 0.8])
%       - EdgeColor   : RGB triplet for polygon edges (default = [0.4 0.4 0.4])
%       - EdgeAlpha   : edge transparency, scalar in [0,1] (default = 0.4)
%       - FaceAlpha   : face transparency, scalar in [0,1] (default = 1)
%       - AxesOpt     : axes handle; if provided, draw on this axes instead of creating a new figure
%       - ValueLabel  : string used as colorbar label (default = 'Value')
% 
% OUTPUT
%   none - this is a renderer that draws directly to a figure or provided axes.
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

if ~isnumeric(Cmap) || ~ismatrix(Cmap) || isempty(Cmap) || size(Cmap,2) ~= 3 || ...
        ~isreal(Cmap) || any(~isfinite(Cmap(:))) || any(Cmap(:) < 0 | Cmap(:) > 1)
    error('options.Colormap must be a nonempty finite M-by-3 numeric matrix with values in [0,1].');
end
if ~isnumeric(NaNColor) || numel(NaNColor) ~= 3 || ~isreal(NaNColor) || ...
        any(~isfinite(NaNColor(:))) || any(NaNColor(:) < 0 | NaNColor(:) > 1)
    error('options.NaNColor must be a finite RGB triplet with values in [0,1].');
end
NaNColor = reshape(NaNColor,1,3);
if ~isnumeric(EdgeColor) || numel(EdgeColor) ~= 3 || ~isreal(EdgeColor) || ...
        any(~isfinite(EdgeColor(:))) || any(EdgeColor(:) < 0 | EdgeColor(:) > 1)
    error('options.EdgeColor must be a finite RGB triplet with values in [0,1].');
end
EdgeColor = reshape(EdgeColor,1,3);
if ~isnumeric(EdgeAlpha) || ~isscalar(EdgeAlpha) || ~isfinite(EdgeAlpha) || ...
        ~isreal(EdgeAlpha) || EdgeAlpha < 0 || EdgeAlpha > 1
    error('options.EdgeAlpha must be a finite scalar in [0,1].');
end
if ~isnumeric(FaceAlpha) || ~isscalar(FaceAlpha) || ~isfinite(FaceAlpha) || ...
        ~isreal(FaceAlpha) || FaceAlpha < 0 || FaceAlpha > 1
    error('options.FaceAlpha must be a finite scalar in [0,1].');
end

%% Step 3: figure/axes setup
if isempty(AxesOpt)
    figure('Color',FigureColor);
    ax = axes;
else
    ax = AxesOpt; axes(ax);
end
holdState = ishold(ax);
hold(ax,'on');

if isempty(CLim)
    finiteVals = Value(isfinite(Value));
    if isempty(finiteVals)
        CLim = [0 1];
    else
        CLim = [min(finiteVals) max(finiteVals)];
        if CLim(1)==CLim(2), CLim = CLim + [-0.5 0.5]; end
    end
else
    if ~isnumeric(CLim) || numel(CLim) ~= 2 || ~isreal(CLim) || any(~isfinite(CLim(:)))
        error('options.ValueLimits must be a finite 1-by-2 numeric vector [min max].');
    end
    CLim = CLim(:).';
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
    if ~isfinite(v)
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
axis(ax,'off');

cb = colorbar(ax);
cb.FontSize = 16;
cb.FontName = 'Times New Roman';
cb.Label.String = ValueLabel;  
cb.Label.FontSize = 16;  
cb.Label.FontName = 'Times New Roman';

if ~holdState, hold(ax,'off'); end
end

function out = gopt(s,name,default)
if isfield(s,name) && ~isempty(s.(name)), out = s.(name); else, out = default; end
end

function cmap2 = resampleColormap(cmap, targetN)
m = size(cmap,1);
if m == 1
    cmap2 = repmat(cmap,targetN,1);
    return;
end
xi = linspace(1,m,m);
xo = linspace(1,m,targetN);
cmap2 = [interp1(xi,cmap(:,1),xo,'linear').', ...
    interp1(xi,cmap(:,2),xo,'linear').', ...
    interp1(xi,cmap(:,3),xo,'linear').'];
cmap2 = max(0,min(1,cmap2));
end
