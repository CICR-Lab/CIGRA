function com_set = covered_component_set_3D(rs, center, radius)
% 输入:
%   rs - 系统结构体，包含 Node(i).Coord3D, Edge(e).FromNodeID, ToNodeID, Dist, CenterCoord3D
%   center - 攻击中心坐标 (1×3)
%   radius - 攻击半径
% 输出:
%   com_set - 所有在攻击范围内的节点/边组成的组件集合（索引：1~N为节点，N+1~N+E为边）

R = 6378.137;
N = length(rs.Node);
E = length(rs.Edge);

% 计算节点与中心点之间的球面距离
node_coords = reshape([rs.Node.Coord3D], 3, [])'; % N×3
ndist_vec = sqrt(sum((node_coords - center).^2, 2));     % 欧几里得距离
ndist = 2 * R * asin(ndist_vec / (2 * R));               % 转为球面距离

% 找到被覆盖的节点索引
com_set = find(ndist <= radius + 1e-5);

% 找到距离中心最近的节点索引
[~, mi] = min(ndist);
mn = min(ndist);

% 准备边起终点对应的节点距离
from_ids = [rs.Edge.FromNodeID];
to_ids   = [rs.Edge.ToNodeID];

from_dist = ndist(from_ids);
to_dist   = ndist(to_ids);

% 距离中值计算
medist1 = min(from_dist, to_dist);  % 任一端在圆内
medist2 = min((from_dist + to_dist - [rs.Edge.Length]') / 2, medist1); % 保守估计
medist3 = max(rs.Dist(N+1:N+E, mi) - mn, medist2);  % 修正值

% 三种方式取最小距离
edist = [(1:E)', medist3, medist1];

% 边直接被包含：中心距离小于 radius
com_set = [com_set; edist(edist(:,3) <= radius + 1e-5, 1) + N];

% 临界边：端点在圆内但中点在边界
temp = edist(edist(:,3) > radius + 1e-5 & edist(:,2) <= radius + 1e-5, 1);

ecom_set = zeros(N, 1); etag = 0;
for t = 1:length(temp)
    edge_idx = temp(t);
    edge_center = rs.Edge(edge_idx).CenterCoord3D;
    edge_dist = norm(edge_center - center);  % 欧氏距离
    if edge_dist <= radius + 1e-5
        etag = etag + 1;
        ecom_set(etag) = N + edge_idx;
    end
end

com_set = [com_set; ecom_set(1:etag)];
end
