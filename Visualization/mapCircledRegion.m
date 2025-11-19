function mapCircledRegion(Boundary, CIS, Mode, layoutOptions, Region)
% INTRODUCTION:
%   mapCircledRegion overlays one or multiple user-specified regions on top
%   of a CIS layout.
%   It first calls mapSingleCISLayout(...) to render the network, then draws
%   one *or more* regions of the following shapes:
%     - 'circle'    : (Center, Radius, Units)
%     - 'polygon'   : (PolygonX, PolygonY)
%     - 'rectangle' : (Center, Width, Height, Units, OrientationDeg)
%
% INPUTS:
%   Boundary, CIS, Mode, layoutOptions : forwarded to mapSingleCISLayout(...)
%
%   Region : either a scalar struct or a 1¡ÁN struct array, each element with
%            fields:
%     Common:
%       - Shape           : 'circle' | 'polygon' | 'rectangle'   (required)
%       - FaceColor       : [r g b], default [1 0 0]
%       - FaceAlpha       : 0~1,     default 0.2
%       - EdgeColor       : [r g b], default [0.2 0.2 0.2]
%       - EdgeAlpha       : 0~1,     default 1.0
%       - EdgeWidth       : numeric, default 1.2
%       - ShowLabel       : true/false, default false
%       - Label           : text, default 'Region'
%       - LabelFontSize   : default 11
%       - LabelColor      : default [0.1 0.1 0.1]
%       - LabelOffset_m   : label offset in meters (north), default 0
%     circle:
%       - Center          : [lon lat]     (required)
%       - Radius          : scalar > 0    (required)
%       - Units           : 'm' | 'km'    (default 'km')
%       - NPoints         : curve smoothness, default 180
%     polygon:
%       - PolygonX, PolygonY : same-length lon/lat vectors (required)
%     rectangle:
%       - Center          : [lon lat]     (required)
%       - Width, Height   : scalars > 0   (required)
%       - Units           : 'm' | 'km'    (default 'km')
%       - OrientationDeg  : rotation w.r.t. east axis, degrees, default 0
%
% OUTPUTS:
%   None (renders into current figure).

%% Step 1: required args & basic checks
if nargin < 5
    error('Region must be a struct or struct array.');
end
if nargin < 4 || isempty(layoutOptions), layoutOptions = struct(); end
if isempty (Region)
    warning ('Empty Region.');
end
%% Step 2: draw base layout only once
mapSingleCISLayout(Boundary, CIS, Mode, layoutOptions);
ax = gca;
hold(ax, 'on');

%% Step 3: loop over regions (support multi-region)
nRegion = numel(Region);
for k = 1:nRegion
    Rk = Region(k);
    draw_single_region(ax, Rk);
end

axis(ax,'equal');
axis(ax,'tight');
hold(ax,'off');

end % main function


% ================= Single-region drawer =================
function draw_single_region(ax, Region)
    % Read common style
    shape          = lower(string(gopt(Region,'Shape','circle')));
    FaceColor      = gopt(Region,'FaceColor',[1 0 0]);
    FaceAlpha      = gopt(Region,'FaceAlpha',0.2);
    EdgeColor      = gopt(Region,'EdgeColor',[0.2 0.2 0.2]);
    EdgeAlpha      = gopt(Region,'EdgeAlpha',1.0);
    EdgeWidth      = gopt(Region,'EdgeWidth',1.2);
    ShowLabel      = gopt(Region,'ShowLabel',false);
    Label          = string(gopt(Region,'Label','Region'));
    LabelFontSize  = gopt(Region,'LabelFontSize',11);
    LabelColor     = gopt(Region,'LabelColor',[0.1 0.1 0.1]);
    LabelOffset_m  = gopt(Region,'LabelOffset_m',0);
    NPoints        = gopt(Region,'NPoints',180); % used by circle

    switch shape
        case "polygon"
            Px = gopt(Region,'PolygonX',[]);
            Py = gopt(Region,'PolygonY',[]);
            if isempty(Px) || isempty(Py) || numel(Px)~=numel(Py)
                error('polygon requires Region.PolygonX/PolygonY of equal length.');
            end
            patch('XData',Px,'YData',Py, ...
                'FaceColor',FaceColor,'FaceAlpha',FaceAlpha, ...
                'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',EdgeWidth, ...
                'Parent',ax);
            ctr = [mean(Px(:),'omitnan'), mean(Py(:),'omitnan')];

        case "circle"
            Center = gopt(Region,'Center',[]);
            if isempty(Center) || numel(Center)~=2
                error('circle requires Region.Center = [lon lat].');
            end
            lon0 = Center(1); lat0 = Center(2);
            Units = lower(string(gopt(Region,'Units','km')));
            Rm    = unit_to_m(gopt(Region,'Radius',[]), Units);
            if isempty(Rm) || ~(Rm>0), error('circle requires positive Radius.'); end

            [dpm_lon, dpm_lat] = deg_per_meter_at_lat(lat0);
            [Px,Py] = circle_ll(lon0, lat0, Rm, NPoints, dpm_lon, dpm_lat);
            patch('XData',Px,'YData',Py, ...
                'FaceColor',FaceColor,'FaceAlpha',FaceAlpha, ...
                'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',EdgeWidth, ...
                'Parent',ax);
            ctr = [lon0, lat0];

        case "rectangle"
            Center = gopt(Region,'Center',[]);
            if isempty(Center) || numel(Center)~=2
                error('rectangle requires Region.Center = [lon lat].');
            end
            lon0 = Center(1); lat0 = Center(2);
            Units = lower(string(gopt(Region,'Units','km')));
            Wm    = unit_to_m(gopt(Region,'Width',[]),  Units);
            Hm    = unit_to_m(gopt(Region,'Height',[]), Units);
            ang   = gopt(Region,'OrientationDeg',0);
            if isempty(Wm) || isempty(Hm) || ~(Wm>0 && Hm>0)
                error('rectangle requires positive Width/Height.');
            end

            [dpm_lon, dpm_lat] = deg_per_meter_at_lat(lat0);
            [Px,Py] = rectangle_ll(lon0, lat0, Wm, Hm, ang, dpm_lon, dpm_lat);
            patch('XData',Px,'YData',Py, ...
                'FaceColor',FaceColor,'FaceAlpha',FaceAlpha, ...
                'EdgeColor',EdgeColor,'EdgeAlpha',EdgeAlpha,'LineWidth',EdgeWidth, ...
                'Parent',ax);
            ctr = [lon0, lat0];

        otherwise
            error('Unsupported Region.Shape = "%s". Use circle|polygon|rectangle.', shape);
    end

    % Label (optional, per-region)
    if ShowLabel
        lonL = ctr(1); latL = ctr(2);
        if LabelOffset_m ~= 0
            [dpm_lon, dpm_lat] = deg_per_meter_at_lat(latL);
            [dlon, dlat] = en_m_to_dlonlat(0, LabelOffset_m, dpm_lon, dpm_lat); % northward offset
            lonL = lonL + dlon; latL = latL + dlat;
        end
        text(lonL, latL, Label, 'Parent',ax, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',LabelFontSize, 'Color',LabelColor, 'FontWeight','bold');
    end
end


% ================= Helpers =================
function v = gopt(s, name, default)
    if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end

function meters = unit_to_m(v, units)
    if isempty(v), meters = []; return; end
    switch string(units)
        case "m",  meters = v;
        case "km", meters = v * 1000;
        otherwise, error('Units must be "m" or "km".');
    end
end

function [deg_per_m_lon, deg_per_m_lat] = deg_per_meter_at_lat(lat)
% Approximate degree-per-meter at given latitude (WGS84 mean values)
    meters_per_deg_lat = 110540;              % ~110.54 km per degree latitude
    meters_per_deg_lon = 111320 * cosd(lat);  % shrinks with latitude
    deg_per_m_lat = 1 / meters_per_deg_lat;
    deg_per_m_lon = 1 / max(1, meters_per_deg_lon); % guard near poles
end

function [dlon, dlat] = en_m_to_dlonlat(east_m, north_m, dpm_lon, dpm_lat)
    dlon = east_m  * dpm_lon;
    dlat = north_m * dpm_lat;
end

function [Px,Py] = circle_ll(lon0, lat0, R_m, N, dpm_lon, dpm_lat)
    t = linspace(0, 2*pi, max(12,N));
    ex = R_m * cos(t); ey = R_m * sin(t);
    [dlon, dlat] = en_m_to_dlonlat(ex, ey, dpm_lon, dpm_lat);
    Px = lon0 + dlon; Py = lat0 + dlat;
end

function [Px,Py] = rectangle_ll(lon0, lat0, W_m, H_m, ang_deg, dpm_lon, dpm_lat)
    hw = W_m/2; hh = H_m/2;
    C = [ -hw, -hh;
           hw, -hh;
           hw,  hh;
          -hw,  hh;
          -hw, -hh ];
    ca = cosd(ang_deg); sa = sind(ang_deg);
    R = [ca -sa; sa ca];
    EN = (R * C.').';
    [dlon, dlat] = en_m_to_dlonlat(EN(:,1), EN(:,2), dpm_lon, dpm_lat);
    Px = lon0 + dlon; Py = lat0 + dlat;
end
