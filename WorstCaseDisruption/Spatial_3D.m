function CIS = Spatial_3D(CIS)
    % Construct 3D Spatial Information
    R = 6378.137; % Earth radius (km)
    N = length(CIS.Node);
    E = length(CIS.Edge);
    
    % Calculate 3D coordinates for nodes
    for m = 1:N
        lon = CIS.Node(m).Longitude;
        lat = CIS.Node(m).Latitude;
        x = R * cosd(lat) * sind(lon);
        y = R * cosd(lat) * cosd(lon);
        z = R * sind(lat);
        CIS.Node(m).Coord3D = [x, y, z];
    end
    
    % Calculate center points for edges
    for e = 1:E
        id1 = CIS.Edge(e).FromNodeID;
        id2 = CIS.Edge(e).ToNodeID;
        coord1 = CIS.Node(id1).Coord3D;
        coord2 = CIS.Node(id2).Coord3D;
        center = (coord1 + coord2) / 2;
        CIS.Edge(e).CenterCoord3D = center;
    end
    
    % Build distance matrix
    dist = zeros(N + E, N + E);
    for i = 1:N + E
        for j = (i + 1):N + E
            % Get coordinates for component i
            if i <= N
                coord_i = CIS.Node(i).Coord3D;
            else
                coord_i = CIS.Edge(i - N).CenterCoord3D;
            end
            
            % Get coordinates for component j
            if j <= N
                coord_j = CIS.Node(j).Coord3D;
            else
                coord_j = CIS.Edge(j - N).CenterCoord3D;
            end
            
            % Calculate Euclidean distance
            dist(i, j) = norm(coord_i - coord_j);
            dist(j, i) = dist(i, j);
        end
    end
    
    % Store distance information
    CIS.Dist = dist;
    for e = 1:E
        id1 = CIS.Edge(e).FromNodeID;
        id2 = CIS.Edge(e).ToNodeID;
        CIS.Edge(e).Dist = dist(id1, id2);
    end
end

