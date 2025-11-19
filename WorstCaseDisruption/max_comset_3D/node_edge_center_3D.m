function [attack_center, initial_coms] = node_edge_center_3D(rs, nodeA, edgeB, radius)
% 输入：
%   rs - 系统结构体（包含 Node(i).Coord3D, Edge(e).FromNodeID/ToNodeID, Edge(e).Length, Edge(e).CenterCoord3D）
%   nodeA - 节点索引
%   edgeB - 边索引
%   radius - 攻击半径（单位：km）
% 输出：
%   attack_center - 中心坐标集合（Nx3）
%   initial_coms - 被初始覆盖的组件编号（节点编号 + N + 边编号）

R = 6378.137;
N = length(rs.Node);

% 如果 nodeA 距离所有其他节点都大于 2r，说明它是孤立点
if sum(rs.Dist(nodeA, 1:N) <= 2 * radius + 1e-5) == 1
    attack_center = rs.Node(nodeA).Coord3D;
    initial_coms = nodeA;
    return;
end

% 如果 nodeA 到 edgeB 的中心距离 > 2r，返回空
if rs.Dist(nodeA, N + edgeB) > 2 * radius + 1e-5
    attack_center = [];
    initial_coms = [];
    return;
end

% 如果 nodeA 是 edgeB 的端点之一，跳过（避免重复）
from_id = rs.Edge(edgeB).FromNodeID;
to_id = rs.Edge(edgeB).ToNodeID;

if ismember(nodeA, [from_id, to_id])
    attack_center = [];
    initial_coms = [];
    return;
end

% 如果 edgeB 是一条极短边（近似为点），转换为两节点处理
if rs.Dist(from_id, to_id) <= 1e-5
    [attack_center, initial_coms] = two_nodes_center_3D(rs, nodeA, from_id, radius);
    return;
end

% 一般情况：通过几何推导计算潜在交点
coord1 = rs.Node(from_id).Coord3D;
coord2 = rs.Node(to_id).Coord3D;

g1 = coord1(2) * coord2(3) - coord2(2) * coord1(3);
g2 = coord2(1) * coord1(3) - coord1(1) * coord2(3);
g3 = coord1(1) * coord2(2) - coord2(1) * coord1(2);
k1 = R * sin(radius / R) * sqrt(1 / (g1^2 + g2^2 + g3^2)) * (g1^2 + g2^2 + g3^2);

D = 2 * R^2 - (2 * R * sin(radius / (2 * R)))^2;
P = 2 * rs.Node(nodeA).Coord3D;

% 解两个球面与大圆交集问题
temp1 = solve_two_curves_equations_3D([g1, g2, g3, k1], [P, D]);
temp2 = solve_two_curves_equations_3D([g1, g2, g3, -k1], [P, D]);
temp = [temp1; temp2];

attack_center = [];
ai = 0;
for t = 1:size(temp, 1)
    dist = norm(rs.Edge(edgeB).CenterCoord3D - temp(t, :));  % 欧几里得距离判断
    if dist <= radius + 1e-5
        ai = ai + 1;
        attack_center(ai, :) = temp(t, :);
    end
end

% 输出结果
if ai > 0
    initial_coms = [nodeA; N + edgeB];
else
    attack_center = [];
    initial_coms = [];
end
end
