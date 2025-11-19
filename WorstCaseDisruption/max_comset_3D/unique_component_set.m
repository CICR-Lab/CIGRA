function ncom_set = unique_component_set(com_set, new_com_set)
% 输入：
%   com_set: 已有的组件集合，维度为 Ns × (3 + Nc)，每行包括攻击中心(1×3)和组件编号集合(1×Nc)
%   new_com_set: 新的候选攻击中心及组件集合 (1 × (3 + Nci))
% 输出：
%   ncom_set: 更新后的 com_set，避免重复组件组合

ntag = 0;

if ~isempty(com_set)
    ncom_set = com_set;
    [N, Nc] = size(com_set);        % N = 当前已有集合数, Nc = 当前最大维度
    Nb = length(new_com_set);       % Nb = 新集合维度

    % 提取组件编号（去掉前3个中心坐标）
    coms = com_set(:, 4:Nc);        % 原集合中每一行的组件编号子集
    Bs = new_com_set(4:Nb);         % 新集合中的组件编号（向量）

    % 在已有集合中查找等价集合
    temp = sum(ismember(coms, Bs), 2);           % 每一行中，与 Bs 重合的数量
    num_non_zero = sum(coms > 0, 2);             % 每一行有效组件个数（非零）
    mi = find(temp == num_non_zero);             % 找到与 Bs 完全重合的行索引

    if ~isempty(mi)
        % 如果存在重复项，则更新第一项，删除其余项
        ntag = 1;
        ncom_set(mi(1), 1:Nb) = new_com_set;
        if length(mi) > 1
            ncom_set(mi(2:end), :) = [];
        end
    else
        % 如果部分重合（组件数一致），也认为是重复
        if any(temp == Nb - 3)
            ntag = 1;
        end
    end

    % 如果是新组件集合，添加进 ncom_set 中
    if ntag == 0
        if Nb <= Nc
            temp = zeros(1, Nc);
            temp(1:Nb) = new_com_set;
            ncom_set = [ncom_set; temp];
        else
            % 补零扩展已有 com_set 到 Nb 列
            ncom_set = [ncom_set, zeros(N, Nb - Nc); new_com_set];
        end
    end

else
    % 初始为空，则直接赋值
    ncom_set = new_com_set;
end
end
