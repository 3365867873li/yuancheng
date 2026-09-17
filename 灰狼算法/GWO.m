function [Alpha_score, Alpha_pos, Convergence_curve] = GWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj)
% GWO 灰狼优化算法主函数
% 输入：SearchAgents_no 种群规模(狼群数量)；Max_iter 最大迭代次数
%       lb/ub 变量上下界；dim 变量维度；fobj 目标函数句柄
% 输出：Alpha_score 最优适应度；Alpha_pos 最优解位置；Convergence_curve 收敛曲线

% 若上下界是标量，则扩展为与维度等长的向量
if numel(ub) == 1, ub = ub .* ones(1, dim); end
if numel(lb) == 1, lb = lb .* ones(1, dim); end

% 初始化 α(头狼)、β(第二狼)、δ(第三狼) 三只领导狼的记录
Alpha_pos = zeros(1, dim); Alpha_score = inf;   % α 狼：当前全局最优
Beta_pos  = zeros(1, dim); Beta_score  = inf;   % β 狼：当前次优
Delta_pos = zeros(1, dim); Delta_score = inf;   % δ 狼：当前第三优

% 第 1 步：随机初始化整个狼群的位置
Positions = initialization(SearchAgents_no, dim, ub, lb);

% 预分配收敛曲线数组，用于记录每次迭代的最优适应度
Convergence_curve = zeros(1, Max_iter);

% ===================== 主迭代循环 =====================
for l = 1:Max_iter

    % 第 2 步：评估每只狼的适应度，并更新 α、β、δ 三只领导狼
    for i = 1:size(Positions, 1)
        % 越界处理：把跑出边界的狼拉回搜索空间
        Positions(i,:) = min(max(Positions(i,:), lb), ub);
        % 计算该狼的适应度
        fitness = fobj(Positions(i,:));

        % 按适应度大小依次更新 α、β、δ
        if fitness < Alpha_score                          % 优于 α
            Alpha_score = fitness; Alpha_pos = Positions(i,:);
        elseif fitness < Beta_score                       % 介于 α 与 β 之间
            Beta_score = fitness; Beta_pos = Positions(i,:);
        elseif fitness < Delta_score                      % 介于 β 与 δ 之间
            Delta_score = fitness; Delta_pos = Positions(i,:);
        end
    end

    % 第 3 步：计算控制参数 a，从 2 线性递减到 0
    % a 控制探索(exploration)到开发(exploitation)的过渡
    a = 2 - l * (2 / Max_iter);

    % 第 4 步：根据 α、β、δ 的位置更新每只狼(ω 狼)的位置
    for i = 1:size(Positions, 1)
        for j = 1:size(Positions, 2)

            % ---------- 受 α 狼引导 ----------
            r1 = rand(); r2 = rand();                 % 随机数 r1,r2 ∈ [0,1]
            A1 = 2*a*r1 - a;                          % 系数 A1，|A|>1 探索，|A|<1 开发
            C1 = 2*r2;                                % 系数 C1，提供随机权重
            D_alpha = abs(C1*Alpha_pos(j) - Positions(i,j)); % 与 α 的距离
            X1 = Alpha_pos(j) - A1*D_alpha;           % 由 α 引导的候选位置

            % ---------- 受 β 狼引导 ----------
            r1 = rand(); r2 = rand();
            A2 = 2*a*r1 - a;
            C2 = 2*r2;
            D_beta = abs(C2*Beta_pos(j) - Positions(i,j)); % 与 β 的距离
            X2 = Beta_pos(j) - A2*D_beta;             % 由 β 引导的候选位置

            % ---------- 受 δ 狼引导 ----------
            r1 = rand(); r2 = rand();
            A3 = 2*a*r1 - a;
            C3 = 2*r2;
            D_delta = abs(C3*Delta_pos(j) - Positions(i,j)); % 与 δ 的距离
            X3 = Delta_pos(j) - A3*D_delta;           % 由 δ 引导的候选位置

            % 取三个候选位置的算术平均作为该维度的新位置
            Positions(i,j) = (X1 + X2 + X3) / 3;
        end
    end

    % 记录本次迭代的全局最优适应度
    Convergence_curve(l) = Alpha_score;
end
end