function [attack_center, initial_coms] = two_nodes_center_3D(rs, nodeA, nodeB, radius)
% Inputs:
%   rs - System structure containing Node(i).Coord3D and the Dist matrix
%   nodeA, nodeB - Node indices
%   radius - Attack radius
% Outputs:
%   attack_center - Attack center coordinates (1-by-3)
%   initial_coms - Indices of the initially covered nodes

R = 6378.137;  % Earth radius
N = length(rs.Node);

% Check whether nodeA is the only node within twice the attack radius.
if sum(rs.Dist(nodeA, 1:N) <= 2 * radius + 1e-5) == 1
    attack_center = rs.Node(nodeA).Coord3D;
    initial_coms = nodeA;

else
    if rs.Dist(nodeA, nodeB) <= 2 * radius + 1e-5
        if rs.Dist(nodeA, nodeB) <= 1e-5
            attack_center = rs.Node(nodeA).Coord3D;
        else
            D = 2 * R^2 - (2 * R * sin(radius / (2 * R)))^2;
            P1 = 2 * rs.Node(nodeA).Coord3D;
            P2 = 2 * rs.Node(nodeB).Coord3D;
            attack_center = solve_two_curves_equations_3D([P1, D], [P2, D]);
        end
        initial_coms = [nodeA; nodeB];
    else
        attack_center = [];
        initial_coms = [];
    end
end
end
