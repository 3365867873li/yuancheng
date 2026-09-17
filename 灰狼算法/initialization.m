function Positions = initialization(SearchAgents_no, dim, ub, lb)
% 随机初始化狼群位置
% 输入：SearchAgents_no 种群规模；dim 维度；ub/lb 上下界
% 输出：Positions 初始位置矩阵(行=狼，列=维度)

if numel(ub) == 1
    % 边界为标量：整块矩阵一次性用均匀随机生成
    Positions = rand(SearchAgents_no, dim) .* (ub - lb) + lb;
else
    % 边界为向量：逐维度在各自 [lb(i), ub(i)] 内均匀采样
    Positions = zeros(SearchAgents_no, dim);
    for i = 1:dim
        Positions(:, i) = rand(SearchAgents_no, 1) .* (ub(i) - lb(i)) + lb(i);
    end
end
end