function [attack_center, initial_coms] = two_nodes_center_3D(rs, nodeA, nodeB, radius)
% 输入：
%   rs - 系统结构体，包含 Node(i).Coord3D 和 Dist 矩阵
%   nodeA, nodeB - 节点编号（索引）
%   radius - 攻击半径
% 输出：
%   attack_center - (1×3) 三维中心坐标
%   initial_coms - 被初步覆盖的组件索引集合（节点编号）

R = 6378.137;  % 地球半径
N = length(rs.Node);

% 如果 nodeA 只有它自己在攻击半径内
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
