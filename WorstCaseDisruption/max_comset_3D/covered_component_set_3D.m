function com_set = covered_component_set_3D(rs, center, radius)
% Inputs:
%   rs - System structure containing node coordinates, edge data, and Dist
%   center - Attack center coordinates (1-by-3)
%   radius - Attack radius
% Output:
%   com_set - Component indices within the attack range; nodes use 1 to N,
%             and edges use N + 1 to N + E

R = 6378.137;
N = length(rs.Node);
E = length(rs.Edge);

% Compute spherical distances from the attack center to all nodes.
node_coords = reshape([rs.Node.Coord3D], 3, [])'; % N-by-3
ndist_vec = sqrt(sum((node_coords - center).^2, 2));     % Euclidean chord distance
ndist = 2 * R * asin(ndist_vec / (2 * R));               % Convert to spherical distance

% Find the covered node indices.
com_set = find(ndist <= radius + 1e-5);

% Find the node nearest to the attack center.
[~, mi] = min(ndist);
mn = min(ndist);

% Get the center-to-node distances at each edge endpoint.
from_ids = [rs.Edge.FromNodeID];
to_ids   = [rs.Edge.ToNodeID];

from_dist = ndist(from_ids);
to_dist   = ndist(to_ids);

% Compute the intermediate edge-distance estimate.
medist1 = min(from_dist, to_dist);  % Distance to the nearer endpoint
medist2 = min((from_dist + to_dist - [rs.Edge.Length]') / 2, medist1); % Conservative estimate
medist3 = max(rs.Dist(N+1:N+E, mi) - mn, medist2);  % Corrected estimate

% Take the minimum of the three edge-distance estimates.
edist = [(1:E)', medist3, medist1];

% Include edges whose estimated distance is within the attack radius.
com_set = [com_set; edist(edist(:,3) <= radius + 1e-5, 1) + N];

% Check the remaining boundary candidates using their edge centers.
temp = edist(edist(:,3) > radius + 1e-5 & edist(:,2) <= radius + 1e-5, 1);

ecom_set = zeros(N, 1); etag = 0;
for t = 1:length(temp)
    edge_idx = temp(t);
    edge_center = rs.Edge(edge_idx).CenterCoord3D;
    edge_dist = norm(edge_center - center);  % Euclidean distance
    if edge_dist <= radius + 1e-5
        etag = etag + 1;
        ecom_set(etag) = N + edge_idx;
    end
end

com_set = [com_set; ecom_set(1:etag)];
end
