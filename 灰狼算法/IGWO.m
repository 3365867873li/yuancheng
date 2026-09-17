function [Alpha_score, Alpha_pos, Convergence_curve, Trigger_log, Alpha_pos_hist] = IGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, cfg)
% ============================================================
% 融合免疫克隆变异机制的改进灰狼算法（论文 2.2 节）
% 在标准 GWO 的位置更新之后，增加"免疫克隆变异"算子：
%   步骤1 计算适应度，更新 α/β/δ
%   步骤2 位置更新（标准灰狼狩猎行为）
%   步骤3 变异触发判断；若触发则对 α/β/δ 克隆→概率变异→精英保留
%   步骤4 参数更新（a 线性衰减）
%   步骤5 结束条件判断
%
% 输入：cfg.immu 提供克隆变异参数（见 mountain_config.m）
% 输出：Trigger_log 记录每代是否触发变异，便于分析
% ============================================================

if numel(ub) == 1, ub = ub .* ones(1, dim); end
if numel(lb) == 1, lb = lb .* ones(1, dim); end

im    = cfg.immu;          % 免疫算子参数
range = ub - lb;           % 各维变量范围，用于缩放变异步长

% 初始化 α、β、δ 三只领导狼
Alpha_pos = zeros(1, dim); Alpha_score = inf;
Beta_pos  = zeros(1, dim); Beta_score  = inf;
Delta_pos = zeros(1, dim); Delta_score = inf;

Positions = initialization(SearchAgents_no, dim, ub, lb);   % 种群初始化
Convergence_curve = zeros(1, Max_iter);                     % 收敛曲线
f_alpha_hist = nan(1, Max_iter);                            % α 狼历史适应度
Trigger_log  = false(1, Max_iter);                          % 触发记录
Alpha_pos_hist = zeros(Max_iter, dim);                      % 逐代最优位置（可选输出）

for g = 1:Max_iter

    %% ---------- 步骤 1：计算适应度，更新 α、β、δ ----------
    for i = 1:size(Positions, 1)
        Positions(i, :) = min(max(Positions(i, :), lb), ub);  % 越界拉回
        fit = fobj(Positions(i, :));
        if fit < Alpha_score
            Alpha_score = fit; Alpha_pos = Positions(i, :);
        elseif fit < Beta_score
            Beta_score = fit;  Beta_pos  = Positions(i, :);
        elseif fit < Delta_score
            Delta_score = fit; Delta_pos = Positions(i, :);
        end
    end
    f_alpha_hist(g) = Alpha_score;

    %% ---------- 步骤 2：位置更新（标准 GWO，式(13)~(15)） ----------
    a = 2 - g * (2 / Max_iter);        % 探索—开发调度因子，线性衰减
    for i = 1:size(Positions, 1)
        for j = 1:dim
            % 受 α 狼引导
            A1 = 2*a*rand() - a;  C1 = 2*rand();
            X1 = Alpha_pos(j) - A1 * abs(C1*Alpha_pos(j) - Positions(i,j));
            % 受 β 狼引导
            A2 = 2*a*rand() - a;  C2 = 2*rand();
            X2 = Beta_pos(j) - A2 * abs(C2*Beta_pos(j) - Positions(i,j));
            % 受 δ 狼引导
            A3 = 2*a*rand() - a;  C3 = 2*rand();
            X3 = Delta_pos(j) - A3 * abs(C3*Delta_pos(j) - Positions(i,j));
            % 三个候选位置的均值作为新位置
            Positions(i, j) = (X1 + X2 + X3) / 3;
        end
    end

    %% ---------- 步骤 3：变异触发判断与免疫克隆变异（论文式(16)） ----------
    % 触发条件：g > gap 且 |f_α(g) - f_α(g-gap)| < ε
    % 含义：前期收敛快、后期停滞时才启动变异，避免无谓地增加计算量
    if g > im.gap && abs(f_alpha_hist(g) - f_alpha_hist(g - im.gap)) < im.eps
        Trigger_log(g) = true;

        leaders = {Alpha_pos, Beta_pos, Delta_pos};       % 待克隆的领导狼
        scores  = [Alpha_score, Beta_score, Delta_score];

        for m = 1:3
            % (1) 克隆：把领导狼复制 Nc 份
            C = repmat(leaders{m}, im.Nc, 1);
            % (2) 概率变异：逐维以概率 pm 施加高斯扰动
            mask = rand(im.Nc, dim) < im.pm;
            C = C + mask .* randn(im.Nc, dim) .* (im.sigma .* range);
            % 边界处理
            C = min(max(C, lb), ub);
            % (3) 评估克隆体，按精英保留策略更新该领导狼
            for k = 1:im.Nc
                fk = fobj(C(k, :));
                if fk < scores(m)
                    scores(m)  = fk;
                    leaders{m} = C(k, :);
                end
            end
        end

        % (4) 从三个成熟抗体群中选出新的 α、β、δ
        [scores, ord] = sort(scores);
        Alpha_score = scores(1); Alpha_pos = leaders{ord(1)};
        Beta_score  = scores(2); Beta_pos  = leaders{ord(2)};
        Delta_score = scores(3); Delta_pos = leaders{ord(3)};
    end

    %% ---------- 步骤 4：记录本次迭代最优 ----------
    Convergence_curve(g) = Alpha_score;
    Alpha_pos_hist(g, :) = Alpha_pos;

end
end
