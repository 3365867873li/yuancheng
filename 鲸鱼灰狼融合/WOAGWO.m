function [Best_score, Best_pos, Convergence_curve] = WOAGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj)
% ============================================================
% 鲸鱼-灰狼融合算法（WOAGWO，Mohammed & Rashid, 2020）
%   论文：A novel hybrid GWO with WOA for global numerical optimization
%         and solving pressure vessel design. Neural Computing and
%         Applications, 32, 14701-14718.
%
% 融合思想（论文 Algorithm 3）：以 WOA 为主框架，把 GWO 的
% α/β/δ 三头狼狩猎机制嵌入 WOA 的开发阶段：
%
%   分支1 p < 0.5（WOA 探索阶段）：
%         |A| < 1   → WOA 包围猎物 Eq.(1)：X = X* − A·|C·X* − X|
%         |A| ≥ 1   → WOA 随机搜索 Eq.(2)：X = X_rand − A·|C·X_rand − X|
%         两个子分支都做"贪婪接受"：仅当新位置适应度更优才更新
%   分支2 p ≥ 0.5（开发阶段）：
%         若 A1、A2、A3 ∈ (−1,1)（GWO 收敛态）：
%           Eq.(14) Dk = |Ck·Xk − X|（k = α, β, δ）
%           Eq.(15) Xk = Xk − Ak·Dk
%           Eq.(16) X = (X1 + X2 + X3) / 3
%         条件不满足则该个体保持原位置（伪代码无 else 分支）
%
% ★ 说明：论文正文描述"开发阶段用 GWO 替换 WOA 螺旋更新，|A|≥1 时
%   仍用螺旋"，但 Algorithm 3 伪代码中开发阶段只有 GWO 式(16) 且
%   无螺旋回退——本实现以 Algorithm 3（正式伪代码）为准。
%
% 参数：a 线性 2→0；A = 2a·r − a；C = 2r；p 均匀随机数
% 输入：SearchAgents_no 种群规模；Max_iter 最大迭代次数；
%       lb/ub 变量下/上界（标量或向量）；dim 维度；fobj 目标函数句柄
% 输出：Best_score 最优适应度；Best_pos 最优解；Convergence_curve 收敛曲线
% ============================================================

if numel(ub) == 1, ub = ub .* ones(1, dim); end
if numel(lb) == 1, lb = lb .* ones(1, dim); end

%% ---------- 初始化种群 ----------
Positions = initialization(SearchAgents_no, dim, ub, lb);  % 随机初始化

% GWO 三头狼：α 最优、β 次优、δ 第三优（每代按适应度更新）
Alpha_pos = zeros(1, dim); Alpha_score = inf;
Beta_pos  = zeros(1, dim); Beta_score  = inf;
Delta_pos = zeros(1, dim); Delta_score = inf;

Convergence_curve = zeros(1, Max_iter);   % 记录每代最优适应度

%% ---------- 主循环 ----------
for t = 1:Max_iter
    a = 2 - t * (2 / Max_iter);   % 收敛系数，从 2 线性递减到 0

    for i = 1:size(Positions, 1)
        % ---- 更新每只个体的 a、A、C、l、p ----
        r1 = rand();  r2 = rand();
        A  = 2*a*r1 - a;   % WOA 系数 A（探索/包围用）
        C  = 2*r2;         % WOA 系数 C
        p  = rand();       % 分支概率：探索/开发各 50%

        % GWO 三头狼各自的 A、C 系数（开发阶段用）
        r1 = rand();  r2 = rand();
        A1 = 2*a*r1 - a;  C1 = 2*r2;   % α 狼系数
        r1 = rand();  r2 = rand();
        A2 = 2*a*r1 - a;  C2 = 2*r2;   % β 狼系数
        r1 = rand();  r2 = rand();
        A3 = 2*a*r1 - a;  C3 = 2*r2;   % δ 狼系数

        f_old = fobj(Positions(i, :));  % 旧位置适应度（贪婪接受用）

        if p < 0.5
            % ================= 分支1：WOA 探索阶段 =================
            Xnew = Positions(i, :);
            if abs(A) < 1
                % ---- 包围猎物 Eq.(1)：向当前最优 X* 靠拢 ----
                for j = 1:dim
                    D = abs(C * Alpha_pos(j) - Positions(i, j));
                    Xnew(j) = Alpha_pos(j) - A * D;
                end
            else
                % ---- 随机搜索 Eq.(2)：随机选一条鲸鱼引导 ----
                rand_idx = randi(SearchAgents_no);   % 随机个体编号
                X_rand   = Positions(rand_idx, :);
                for j = 1:dim
                    D = abs(C * X_rand(j) - Positions(i, j));
                    Xnew(j) = X_rand(j) - A * D;
                end
            end
            % ---- 贪婪接受：仅当新位置适应度更优才更新 ----
            Xnew = min(max(Xnew, lb), ub);            % 越界拉回
            if fobj(Xnew) < f_old
                Positions(i, :) = Xnew;
            end

        else
            % ================= 分支2：开发阶段（GWO 三头狼） =================
            % 式(16) 触发条件：A1、A2、A3 均在 (−1, 1)
            % （伪代码原文为 |Ai| > -1 || |Ai| < 1，因 |Ai| ≥ 0 > -1 恒真，
            %   等价于 |Ai| < 1，即三头狼都处于收敛引导状态）
            if abs(A1) < 1 && abs(A2) < 1 && abs(A3) < 1
                Xnew = Positions(i, :);
                for j = 1:dim
                    % Eq.(14)：与三头狼的距离
                    D1 = abs(C1 * Alpha_pos(j) - Positions(i, j));
                    D2 = abs(C2 * Beta_pos(j)  - Positions(i, j));
                    D3 = abs(C3 * Delta_pos(j) - Positions(i, j));
                    % Eq.(15)：三头狼各自的引导位置
                    X1 = Alpha_pos(j) - A1 * D1;
                    X2 = Beta_pos(j)  - A2 * D2;
                    X3 = Delta_pos(j) - A3 * D3;
                    % Eq.(16)：取三者均值作为新位置
                    Xnew(j) = (X1 + X2 + X3) / 3;
                end
                Positions(i, :) = min(max(Xnew, lb), ub);   % 越界拉回
            end
            % 条件不满足：该个体保持原位置不动（伪代码无 else）
        end
    end

    % ---------- 越界修正、评估适应度、更新 α/β/δ ----------
    for i = 1:size(Positions, 1)
        Positions(i, :) = min(max(Positions(i, :), lb), ub);  % 拉回边界内
        fit = fobj(Positions(i, :));
        if fit < Alpha_score
            Alpha_score = fit;  Alpha_pos = Positions(i, :);
        elseif fit < Beta_score
            Beta_score = fit;   Beta_pos  = Positions(i, :);
        elseif fit < Delta_score
            Delta_score = fit;  Delta_pos = Positions(i, :);
        end
    end

    Convergence_curve(t) = Alpha_score;        % 记录本代最优值
end

Best_score = Alpha_score;
Best_pos   = Alpha_pos;
end
