function zoneInfo = identifyEdgeStrPassingZones(TerminalZone, EdgeStr)
% getZonesLengthFromEdge computes, for each zone, the total length of an edge 
% (given by its turning-point coordinates) that lies within that zone.
%
%   zoneInfo = getZonesLengthFromEdge(TerminalZone, EdgeStr)
%
%   Inputs:
%     TerminalZone - A structure array with fields:
%           TerminalZone(k).X : Vector of longitudes defining the zone boundary 
%                               (with NaN at the end).
%           TerminalZone(k).Y : Vector of latitudes defining the zone boundary
%                               (with NaN at the end).
%           TerminalZone(k).Type : (land use type; can be empty)
%           TerminalZone(k).Population : Population assigned to the zone.
%
%     EdgeStr - A structure with fields:
%           EdgeStr.X : Vector of longitudes along the edge turning points (with NaN at the end)
%           EdgeStr.Y : Vector of latitudes along the edge turning points (with NaN at the end)
%
%   Output:
%     zoneInfo - A structure array where each element has fields:
%           zoneID      : The index of the TerminalZone in which the edge has nonzero coverage.
%           insideLength: Total length (in km) of the portion of the edge that lies inside the zone.
%
%   Method:
%     1. Remove NaN entries from the edge coordinate vectors to form a continuous polyline.
%     2. For each TerminalZone:
%           - Remove NaN values from the zone boundary.
%           - For each consecutive pair of turning points in the edge:
%                 * If both endpoints are inside the polygon, add the full segment length.
%                 * If one endpoint is inside and one outside, use polyxpoly to compute the intersection,
%                   and add the partial segment length inside.
%                 * If both endpoints are outside but the segment crosses the boundary,
%                   use the two intersection points.
%     3. Return the zone IDs (indices) that the edge touches and the total inside length for each.
%
%   Author: [Your Name]
%   Date: [Date]

%% Preprocess EdgeStr: remove NaN values
edgeX = EdgeStr.X;
edgeY = EdgeStr.Y;
validEdge = ~isnan(edgeX) & ~isnan(edgeY);
edgeX = edgeX(validEdge);
edgeY = edgeY(validEdge);

% Conversion factor: 1 degree ~ 111 km (approximate)
zoneInfo = struct('zoneID', {}, 'insideLength', {});

%% Loop over each TerminalZone
for k = 1:length(TerminalZone)
    % Get the polygon for zone k: remove NaN values.
    zoneX = TerminalZone(k).X;
    zoneY = TerminalZone(k).Y;
    validZone = ~isnan(zoneX) & ~isnan(zoneY);
    zoneX = zoneX(validZone);
    zoneY = zoneY(validZone);
    
    % Initialize cumulative inside length.
    insideLength = 0;
    
    % Process each segment of the edge.
    nPts = length(edgeX);
    for i = 1:(nPts-1)
        p1 = [edgeX(i), edgeY(i)];
        p2 = [edgeX(i+1), edgeY(i+1)];
        
        % Check if endpoints are inside the zone.
        in1 = inpolygon(p1(1), p1(2), zoneX, zoneY);
        in2 = inpolygon(p2(1), p2(2), zoneX, zoneY);
        
        if in1 && in2
            % Both endpoints inside: add full segment length.
            segLength = longitude_latitude(p1(1),p1(2),p2(1),p2(2));
            insideLength = insideLength + segLength;
        elseif in1 && ~in2
            % p1 is inside, p2 is outside: find intersection.
            [xi, yi] = polyxpoly([p1(1) p2(1)], [p1(2) p2(2)], zoneX, zoneY);
            if ~isempty(xi)
                pInt = [xi(1), yi(1)];
                segLength = longitude_latitude(p1(1),p1(2),pInt(1),pInt(2));
                insideLength = insideLength + segLength;
            end
        elseif ~in1 && in2
            % p1 is outside, p2 is inside: find intersection.
            [xi, yi] = polyxpoly([p1(1) p2(1)], [p1(2) p2(2)], zoneX, zoneY);
            if ~isempty(xi)
                pInt = [xi(1), yi(1)];
                segLength = longitude_latitude(p2(1),p2(2),pInt(1),pInt(2));
                insideLength = insideLength + segLength;
            end
        else
            % Both endpoints outside. Check if segment crosses the polygon.
            [xi, yi] = polyxpoly([p1(1) p2(1)], [p1(2) p2(2)], zoneX, zoneY);
            if numel(xi) >= 2
                % Use the first two intersections.
                pInt1 = [xi(1), yi(1)];
                pInt2 = [xi(2), yi(2)];
                segLength = longitude_latitude(pInt1(1),pInt1(2),pInt2(1),pInt2(2));
                insideLength = insideLength + segLength;
            end
        end
    end  % end for segments
    
    % If any length was accumulated, then record the zone.
    if insideLength > 0
        zoneInfo(end+1).zoneID = TerminalZone(k).ID;     %#ok<AGROW>
        zoneInfo(end).insideLength = insideLength;
    end
end

end
