function zoneInfo = identifyEdgeStrPassingZones(TerminalZone, EdgeStr)
% identifyEdgeStrPassingZones Identifies zones crossed by an infrastructure edge.
%
%   zoneInfo = identifyEdgeStrPassingZones(TerminalZone, EdgeStr) returns
%   each terminal zone containing a nonzero portion of the edge and the
%   length of that portion in km.
%
%   TerminalZone is a structure array of rectangular zones with fields:
%       ID - zone identifier
%       X  - zone boundary longitudes, with NaN at the end
%       Y  - zone boundary latitudes, with NaN at the end
%
%   EdgeStr is a continuous polyline with fields:
%       X - turning-point longitudes, with NaN at the end
%       Y - turning-point latitudes, with NaN at the end
%
%   zoneInfo is a structure array with fields:
%       zoneID       - value of TerminalZone(k).ID
%       insideLength - total edge length inside that zone, in km
%
%   The function first rejects zones whose bounding boxes do not overlap the
%   edge. It then clips each edge segment against the remaining zone
%   boundaries and accumulates the length inside each zone.

%% Preprocess EdgeStr: remove NaN values
edgeX = EdgeStr.X;
edgeY = EdgeStr.Y;
validEdge = ~isnan(edgeX) & ~isnan(edgeY);
edgeX = edgeX(validEdge);
edgeY = edgeY(validEdge);

zoneInfo = struct('zoneID', {}, 'insideLength', {});
if numel(edgeX) < 2
    return;
end

edgeMinX = min(edgeX);
edgeMaxX = max(edgeX);
edgeMinY = min(edgeY);
edgeMaxY = max(edgeY);

%% Loop over each TerminalZone
for k = 1:length(TerminalZone)
    % Get the polygon for zone k: remove NaN values.
    zoneX = TerminalZone(k).X;
    zoneY = TerminalZone(k).Y;
    validZone = ~isnan(zoneX) & ~isnan(zoneY);
    zoneX = zoneX(validZone);
    zoneY = zoneY(validZone);

    % Skip zones whose bounding boxes do not overlap the edge.
    if isempty(zoneX) || edgeMaxX < min(zoneX) || edgeMinX > max(zoneX) || ...
            edgeMaxY < min(zoneY) || edgeMinY > max(zoneY)
        continue;
    end
    
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
                distanceSquared = (xi - p1(1)).^2 + (yi - p1(2)).^2;
                [~, farthestIdx] = max(distanceSquared);
                pInt = [xi(farthestIdx), yi(farthestIdx)];
                segLength = longitude_latitude(p1(1),p1(2),pInt(1),pInt(2));
                insideLength = insideLength + segLength;
            end
        elseif ~in1 && in2
            % p1 is outside, p2 is inside: find intersection.
            [xi, yi] = polyxpoly([p1(1) p2(1)], [p1(2) p2(2)], zoneX, zoneY);
            if ~isempty(xi)
                distanceSquared = (xi - p2(1)).^2 + (yi - p2(2)).^2;
                [~, farthestIdx] = max(distanceSquared);
                pInt = [xi(farthestIdx), yi(farthestIdx)];
                segLength = longitude_latitude(p2(1),p2(2),pInt(1),pInt(2));
                insideLength = insideLength + segLength;
            end
        else
            % Both endpoints outside. Check if segment crosses the polygon.
            [xi, yi] = polyxpoly([p1(1) p2(1)], [p1(2) p2(2)], zoneX, zoneY);
            if numel(xi) >= 2
                % Use the first and last intersections along the segment.
                segmentDirection = p2 - p1;
                intersectionPosition = (xi - p1(1)) * segmentDirection(1) + ...
                    (yi - p1(2)) * segmentDirection(2);
                [~, firstIdx] = min(intersectionPosition);
                [~, lastIdx] = max(intersectionPosition);
                pInt1 = [xi(firstIdx), yi(firstIdx)];
                pInt2 = [xi(lastIdx), yi(lastIdx)];
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
