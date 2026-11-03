function [attack_center, initial_coms] = two_edges_center_3D(rs, edgeA, edgeB, radius)
% Inputs:
%   rs - System structure containing node coordinates and edge data
%   edgeA, edgeB - Edge indices
%   radius - Attack radius
% Outputs:
%   attack_center - Feasible attack center coordinates (N-by-3)
%   initial_coms - Covered component indices (N + edgeA and N + edgeB)

R = 6378.137;
N = length(rs.Node);

% Check whether the two edge centers are within twice the attack radius.
if rs.Dist(N + edgeA, N + edgeB) <= 2 * radius + 1e-3
    if rs.Dist(N + edgeA, N + edgeB) <= 1e-5
        attack_center = [];
        initial_coms = [];
        return;
    end

    % Initialize the candidate arrays.
    attack_center = zeros(4, 3); ai = 0;

    % Get the endpoint coordinates of edgeA.
    n1 = rs.Edge(edgeA).FromNodeID;
    n2 = rs.Edge(edgeA).ToNodeID;
    coord1 = rs.Node(n1).Coord3D;
    coord2 = rs.Node(n2).Coord3D;

    % Get the endpoint coordinates of edgeB.
    n3 = rs.Edge(edgeB).FromNodeID;
    n4 = rs.Edge(edgeB).ToNodeID;
    coord3 = rs.Node(n3).Coord3D;
    coord4 = rs.Node(n4).Coord3D;

    % Compute the great-circle plane parameters for edgeA.
    g1 = coord1(2) * coord2(3) - coord2(2) * coord1(3);
    g2 = coord2(1) * coord1(3) - coord1(1) * coord2(3);
    g3 = coord1(1) * coord2(2) - coord2(1) * coord1(2);
    k1 = R * sin(radius / R) * sqrt(1 / (g1^2 + g2^2 + g3^2)) * (g1^2 + g2^2 + g3^2);

    % Compute the great-circle plane parameters for edgeB.
    h1 = coord3(2) * coord4(3) - coord4(2) * coord3(3);
    h2 = coord4(1) * coord3(3) - coord3(1) * coord4(3);
    h3 = coord3(1) * coord4(2) - coord4(1) * coord3(2);
    k2 = R * sin(radius / R) * sqrt(1 / (h1^2 + h2^2 + h3^2)) * (h1^2 + h2^2 + h3^2);

    % Solve the four possible pairs of offset great-circle planes.
    temp1 = solve_two_curves_equations_3D([g1, g2, g3, k1],  [h1, h2, h3,  k2]);
    temp2 = solve_two_curves_equations_3D([g1, g2, g3, k1],  [h1, h2, h3, -k2]);
    temp3 = solve_two_curves_equations_3D([g1, g2, g3, -k1], [h1, h2, h3,  k2]);
    temp4 = solve_two_curves_equations_3D([g1, g2, g3, -k1], [h1, h2, h3, -k2]);
    temp = [temp1; temp2; temp3; temp4];

    % Retain centers within the attack radius of both edge centers.
    for t = 1:size(temp, 1)
        center = temp(t, :);
        distA = norm(rs.Edge(edgeA).CenterCoord3D - center);
        distB = norm(rs.Edge(edgeB).CenterCoord3D - center);
        if distA <= radius + 1e-4 && distB <= radius + 1e-4
            ai = ai + 1;
            attack_center(ai, :) = center;
        end
    end

    attack_center = attack_center(1:ai, :);
    initial_coms = [N + edgeA; N + edgeB];

else
    attack_center = [];
    initial_coms = [];
end
end
