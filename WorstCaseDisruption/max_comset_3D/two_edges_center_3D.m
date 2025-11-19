function [attack_center, initial_coms] = two_edges_center_3D(rs, edgeA, edgeB, radius)
% 输入：
%   rs - 系统结构体，包含 Edge(i).FromNodeID/ToNodeID，Node(i).Coord3D，Edge(i).CenterCoord3D
%   edgeA, edgeB - 边的编号（索引）
%   radius - 攻击半径
% 输出：
%   attack_center - 可行攻击中心点集合（Nx3）
%   initial_coms - 被覆盖的组件索引集合（N+edgeA，N+edgeB）

R = 6378.137;
N = length(rs.Node);

% 两边中心距离是否在攻击半径之内
if rs.Dist(N + edgeA, N + edgeB) <= 2 * radius + 1e-3
    if rs.Dist(N + edgeA, N + edgeB) <= 1e-5
        attack_center = [];
        initial_coms = [];
        return;
    end

    % 初始化
    attack_center = zeros(4, 3); ai = 0;

    % 获取边A的两个端点坐标
    n1 = rs.Edge(edgeA).FromNodeID;
    n2 = rs.Edge(edgeA).ToNodeID;
    coord1 = rs.Node(n1).Coord3D;
    coord2 = rs.Node(n2).Coord3D;

    % 获取边B的两个端点坐标
    n3 = rs.Edge(edgeB).FromNodeID;
    n4 = rs.Edge(edgeB).ToNodeID;
    coord3 = rs.Node(n3).Coord3D;
    coord4 = rs.Node(n4).Coord3D;

    % 计算边A的大圆参数
    g1 = coord1(2) * coord2(3) - coord2(2) * coord1(3);
    g2 = coord2(1) * coord1(3) - coord1(1) * coord2(3);
    g3 = coord1(1) * coord2(2) - coord2(1) * coord1(2);
    k1 = R * sin(radius / R) * sqrt(1 / (g1^2 + g2^2 + g3^2)) * (g1^2 + g2^2 + g3^2);

    % 边B的大圆参数
    h1 = coord3(2) * coord4(3) - coord4(2) * coord3(3);
    h2 = coord4(1) * coord3(3) - coord3(1) * coord4(3);
    h3 = coord3(1) * coord4(2) - coord4(1) * coord3(2);
    k2 = R * sin(radius / R) * sqrt(1 / (h1^2 + h2^2 + h3^2)) * (h1^2 + h2^2 + h3^2);

    % 求两大圆交点（共四组可能）
    temp1 = solve_two_curves_equations_3D([g1, g2, g3, k1],  [h1, h2, h3,  k2]);
    temp2 = solve_two_curves_equations_3D([g1, g2, g3, k1],  [h1, h2, h3, -k2]);
    temp3 = solve_two_curves_equations_3D([g1, g2, g3, -k1], [h1, h2, h3,  k2]);
    temp4 = solve_two_curves_equations_3D([g1, g2, g3, -k1], [h1, h2, h3, -k2]);
    temp = [temp1; temp2; temp3; temp4];

    % 筛选在两个边攻击半径范围内的交点
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
