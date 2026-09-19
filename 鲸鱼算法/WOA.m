function [Best_score, Best_pos, Convergence_curve] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj)
% ============================================================
% 鲸鱼优化算法（WOA，Mirjalili & Lewis, 2016）
%   论文：The Whale Optimization Algorithm, Advances in Engineering Software, 95, 51-67
%
% 模拟座头鲸捕食行为的三种策略：
%   策略1 包围猎物   ：D = |C·X* − X|，X(t+1) = X* − A·D
%   策略2 气泡网攻击 ：两种机制按概率 p 切换
%         (a) 收缩包围：a 线性递减使 A ∈ [−a, a] 缩小，p < 0.5 且 |A| < 1
%         (b) 螺旋更新：X(t+1) = D'·e^(b·l)·cos(2πl) + X*，p ≥ 0.5
%   策略3 随机搜索  ：|A| ≥ 1 时随机选择鲸鱼 Xrand 引导，增强全局探索
%
% 内部参数：A = 2a·r − a（|A|<1 开发、≥1 探索），C = 2r，
%          a 从 2 线性递减到 0，b = 1 螺旋形状常数，l ∈ [−1, 1]
%
% 输入：SearchAgents_no 种群规模；Max_iter 最大迭代次数；
%       lb/ub 变量下/上界（标量或向量）；dim 维度；fobj 目标函数句柄
% 输出：Best_score 最优适应度；Best_pos 最优解；Convergence_curve 收敛曲线
% ============================================================

if numel(ub) == 1, ub = ub .* ones(1, dim); end
if numel(lb) == 1, lb = lb .* ones(1, dim); end

%% ---------- 初始化鲸鱼种群 ----------
Positions = initialization(SearchAgents_no, dim, ub, lb);  % 随机初始化

Leader_pos   = zeros(1, dim);   % 领导者（当前最优鲸鱼）位置
Leader_score = inf;             % 领导者适应度（极小化问题，初始为 +∞）

Convergence_curve = zeros(1, Max_iter);   % 记录每代最优适应度

%% ---------- 主循环 ----------
for t = 1:Max_iter
    a = 2 - t * (2 / Max_iter);   % 收敛系数，从 2 线性递减到 0（探索→开发）
    b = 1;                        % 螺旋形状常数（论文取 1）

    for i = 1:size(Positions, 1)
        p  = rand();              % 概率参数，决定收缩包围还是螺旋更新（各 50%）
        r1 = rand();  r2 = rand();
        A  = 2*a*r1 - a;          % 系数 A：|A|≥1 探索，|A|<1 开发
        C  = 2*r2;                % 系数 C
        l  = 2*rand() - 1;        % 螺旋参数 l ∈ [−1, 1]

        for j = 1:dim
            if p < 0.5
                % ---------- 策略1/3：包围猎物或随机搜索 ----------
                if abs(A) >= 1
                    % 探索阶段：随机选一条鲸鱼作为引导，跳出局部最优
                    rand_idx = randi(SearchAgents_no);   % 随机个体编号
                    X_rand   = Positions(rand_idx, :);
                    D        = abs(C * X_rand(j) - Positions(i, j));
                    Positions(i, j) = X_rand(j) - A * D;
                else
                    % 开发阶段：收缩包围，向当前最优鲸鱼靠拢
                    D = abs(C * Leader_pos(j) - Positions(i, j));
                    Positions(i, j) = Leader_pos(j) - A * D;
                end
            else
                % ---------- 策略2(b)：螺旋更新（气泡网攻击） ----------
                Dp = abs(Leader_pos(j) - Positions(i, j));         % 与最优解的距离
                Positions(i, j) = Dp * exp(b*l) * cos(2*pi*l) + Leader_pos(j);
            end
        end
    end

    % ---------- 越界修正、评估适应度、更新领导者 ----------
    for i = 1:size(Positions, 1)
        Positions(i, :) = min(max(Positions(i, :), lb), ub);   % 拉回边界内
        fit = fobj(Positions(i, :));
        if fit < Leader_score
            Leader_score = fit;                 % 更新最优适应度
            Leader_pos   = Positions(i, :);     % 更新最优鲸鱼位置
        end
    end

    Convergence_curve(t) = Leader_score;        % 记录本代最优值
end

Best_score = Leader_score;
Best_pos   = Leader_pos;
end
