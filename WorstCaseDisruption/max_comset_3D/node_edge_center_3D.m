function [attack_center, initial_coms] = node_edge_center_3D(rs, nodeA, edgeB, radius)
% Inputs:
%   rs - System structure containing Node(i).Coord3D and edge data
%   nodeA - Node index
%   edgeB - Edge index
%   radius - Attack radius in km
% Outputs:
%   attack_center - Candidate attack center coordinates (N-by-3)
%   initial_coms - Initially covered component indices (nodeA and N + edgeB)

R = 6378.137;
N = length(rs.Node);

% Treat nodeA as isolated if every other node is farther than 2r.
if sum(rs.Dist(nodeA, 1:N) <= 2 * radius + 1e-5) == 1
    attack_center = rs.Node(nodeA).Coord3D;
    initial_coms = nodeA;
    return;
end

% Return empty if nodeA is farther than 2r from the center of edgeB.
if rs.Dist(nodeA, N + edgeB) > 2 * radius + 1e-5
    attack_center = [];
    initial_coms = [];
    return;
end

% Skip nodeA if it is an endpoint of edgeB to avoid duplication.
from_id = rs.Edge(edgeB).FromNodeID;
to_id = rs.Edge(edgeB).ToNodeID;

if ismember(nodeA, [from_id, to_id])
    attack_center = [];
    initial_coms = [];
    return;
end

% Treat an extremely short edgeB as a point and use the two-node routine.
if rs.Dist(from_id, to_id) <= 1e-5
    [attack_center, initial_coms] = two_nodes_center_3D(rs, nodeA, from_id, radius);
    return;
end

% Compute candidate intersections for the general case.
coord1 = rs.Node(from_id).Coord3D;
coord2 = rs.Node(to_id).Coord3D;

g1 = coord1(2) * coord2(3) - coord2(2) * coord1(3);
g2 = coord2(1) * coord1(3) - coord1(1) * coord2(3);
g3 = coord1(1) * coord2(2) - coord2(1) * coord1(2);
k1 = R * sin(radius / R) * sqrt(1 / (g1^2 + g2^2 + g3^2)) * (g1^2 + g2^2 + g3^2);

D = 2 * R^2 - (2 * R * sin(radius / (2 * R)))^2;
P = 2 * rs.Node(nodeA).Coord3D;

% Solve the intersections of the node-distance plane, the offset great-circle
% planes, and the Earth sphere.
temp1 = solve_two_curves_equations_3D([g1, g2, g3, k1], [P, D]);
temp2 = solve_two_curves_equations_3D([g1, g2, g3, -k1], [P, D]);
temp = [temp1; temp2];

attack_center = [];
ai = 0;
for t = 1:size(temp, 1)
    dist = norm(rs.Edge(edgeB).CenterCoord3D - temp(t, :));  % Euclidean distance check
    if dist <= radius + 1e-5
        ai = ai + 1;
        attack_center(ai, :) = temp(t, :);
    end
end

% Return the candidate centers and their initial component pair.
if ai > 0
    initial_coms = [nodeA; N + edgeB];
else
    attack_center = [];
    initial_coms = [];
end
end
