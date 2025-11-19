function com_set = maximal_component_set_3D(rs, radius)
% For each node-node, node-edge, edge-edge pair, generate attack centers 
% and compute the covered component set within the given radius.

N = length(rs.Node); 
E = length(rs.Edge);
com_set = [];

for i = 1:(N + E)
    for j = (i + 1):(N + E)
        if i <= N
            if j <= N
                [attack_center, initial_coms] = two_nodes_center_3D(rs, i, j, radius);
            else
                [attack_center, initial_coms] = node_edge_center_3D(rs, i, j - N, radius);
            end
        else
            [attack_center, initial_coms] = two_edges_center_3D(rs, i - N, j - N, radius);
        end

        if ~isempty(attack_center)
            for s = 1:size(attack_center, 1)
                cov_com_set = covered_component_set_3D(rs, attack_center(s, :), radius);
                ccom_set = union(cov_com_set, initial_coms);
                new_com_set = [attack_center(s, :) ccom_set'];
                com_set = unique_component_set(com_set, new_com_set);
            end
        end
    end
end
end
