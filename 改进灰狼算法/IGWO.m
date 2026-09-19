function [Alpha_score, Alpha_pos, Convergence_curve, Trigger_log, Alpha_pos_hist] = IGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, cfg)
% ============================================================
% 融合免疫克隆变异机制的改进灰狼算法（论文 2.2 节）+ 两项增强
% 在标准 GWO 的位置更新之后，增加三个改进算子：
%   步骤1 计算适应度，更新 α/β/δ
%   步骤2 位置更新（标准灰狼狩猎行为）
%         · 自适应加权：按适应度给 α/β/δ 分配权重，好狼主导引导
%         · Lévy 飞行：部分狼施加重尾随机跳跃，增强全局探索能力
%   步骤3 免疫克隆变异触发判断（论文式(16)）
%   步骤4 参数更新（a 线性衰减）
%   步骤5 结束条件判断
%
% 输入：cfg.immu 提供克隆变异参数；cfg.levy 提供 Lévy 飞行参数（均可缺省）
% 输出：Trigger_log 记录每代是否触发变异，便于分析
% ============================================================

if numel(ub) == 1, ub = ub .* ones(1, dim); end
if numel(lb) == 1, lb = lb .* ones(1, dim); end

im    = cfg.immu;          % 免疫算子参数
if ~isfield(im, 'minGen'), im.minGen = 30; end   % 兼容未配置该字段的旧调用（论文式(16)：g>30）
range = ub - lb;           % 各维变量范围，用于缩放变异步长

% ---- Lévy 飞行参数（Mantegna 算法生成重尾步长） ----
if ~isfield(cfg, 'levy'), cfg.levy = struct(); end
if ~isfield(cfg.levy, 'p'),     cfg.levy.p      = 0.3;   end  % 每只狼每代施加 Lévy 扰动的概率
if ~isfield(cfg.levy, 'lambda'), cfg.levy.lambda = 1.5;  end  % 重尾指数 λ∈(1,2]，越大跃迁越"重尾"
pLevy   = cfg.levy.p;
lambda  = cfg.levy.lambda;
% Mantegna 算法中的 σ 常数（只与 λ 有关，预先算一次）
sigmaLevy = (gamma(1+lambda) * sin(pi*lambda/2) / ...
             (gamma((1+lambda)/2) * lambda * 2^((lambda-1)/2)))^(1/lambda);

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

    %% ---------- 步骤 2：位置更新（自适应加权 + Lévy 飞行） ----------
    a = 2 - g * (2 / Max_iter);        % 探索—开发调度因子，线性衰减

    % ---- 自适应引导权重：按适应度给 α/β/δ 分配权重（适应度越好权重越大） ----
    % α≤β≤δ 已保证有序；把适应度做 min-max 归一化后整体上移 1，
    % 使三只狼都有非零权重（早期量级差大时 δ 也不会被完全忽略）
    fspan = Delta_score - Alpha_score;                          % 适应度跨度（≥0）
    u     = 1 + (Delta_score - [Alpha_score, Beta_score, Delta_score]) / (fspan + eps);
    w     = u / sum(u);                                         % 归一化：w(1)≥w(2)≥w(3)

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
            % 按自适应权重合成新位置（替代原等权平均 (X1+X2+X3)/3）
            Positions(i, j) = w(1)*X1 + w(2)*X2 + w(3)*X3;
        end

        % ---- Lévy 飞行：以概率 pLevy 对更新后的个体叠加重尾随机步长 ----
        % 步长 ∝ a（前期大、后期小），维持探索—开发平衡，防止陷入局部最优
        if rand < pLevy
            % Mantegna 法：两个高斯随机数之比构成重尾步长
            % 分母加 1e-12 下限，避免 randn 恰好为 0 产生 NaN/Inf 污染种群
            v = abs(randn(1, dim));
            lsteps = randn(1, dim) ./ max(v, 1e-12).^(1/lambda) .* sigmaLevy;
            Positions(i, :) = Positions(i, :) + a .* lsteps .* range;
            Positions(i, :) = min(max(Positions(i, :), lb), ub);   % 越界拉回
        end
    end

    %% ---------- 步骤 3：变异触发判断与免疫克隆变异（论文式(16)） ----------
    % 触发条件：g > 30 且 |f_α(g) - f_α(g-4)| < 0.01（论文式(16)）
    % 含义：前期收敛快、后期停滞时才启动变异，避免无谓地增加计算量
    % 注意：取绝对值——极小化下 f_α 单调不增，f_α(g)-f_α(g-4)≤0 恒成立，
    %       若不取绝对值条件将恒真，失去"仅停滞时变异"的意义
    if g > max(im.minGen, im.gap) && abs(f_alpha_hist(g) - f_alpha_hist(g - im.gap)) < im.eps
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
