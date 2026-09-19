clear; clc; close all;

%% ==================== 用户控制区 ====================
cfg = mountain_config();                    % 载入山地场景与指标配置

% 四项指标权重 [隐蔽性, 机动性, 安全性, 航程]（论文式(12)的 η）
cfg.eta = [0.25, 0.25, 0.25, 0.25];         % <<< 修改这里即可调整偏好

SearchAgents_no = 40;      % 狼群规模
Max_iter        = 500;     % 最大迭代次数
runs            = 10;      % 独立重复运行次数（与论文表1/表2 的 10 次对应）

%% ==================== 构建搜索空间 ====================
n  = cfg.n;   dim = 3 * n;                  % 决策变量维度
xa = cfg.S(1);  xb = cfg.T(1);  dx = (xb - xa) / (n + 1);
lb = zeros(1, dim);   ub = zeros(1, dim);
for i = 1:n
    lb(3*i-2) = xa + (i-1)*dx;   ub(3*i-2) = xa + i*dx;      % x：限制在等分区间
    lb(3*i-1) = cfg.y_range(1);  ub(3*i-1) = cfg.y_range(2); % y
    lb(3*i)   = cfg.z_range(1);  ub(3*i)   = cfg.z_range(2); % z（高度）
end
fobj = @(x) mountain_objective(x, cfg);

fprintf('航路点总数 N_p = %d，决策变量维度 = %d\n', n + 2, dim);
fprintf('权重 [隐蔽 机动 安全 航程] = %s\n\n', mat2str(cfg.eta, 3));

%% ==================== 基线 GWO（无变异） ====================
gwoScore = zeros(runs,1);  gwoTime = zeros(runs,1);      % 最终目标值 / 耗时(ms)
gwoCurve = zeros(runs, Max_iter);                       % 每代惩罚目标值（收敛曲线）
KmatG    = zeros(runs, Max_iter);                       % 逐代加权总指标（式(12)，无罚项）
IndG     = zeros(runs, Max_iter, 4);                    % 逐代四指标 [隐蔽Jr 机动Jm 安全Js 航程Jl]
gwoBest = inf;  gwoX = [];  gwoM = [];

for r = 1:runs
    rng(r);                                 % 固定随机种子，保证两算法配对可比
    t0 = tic;
    [s, pos, conv, posHist] = GWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
    gwoTime(r) = toc(t0) * 1000;            % 记录耗时 / ms（论文表1/表2 单位为 ms）
    gwoScore(r) = s;  gwoCurve(r,:) = conv;
    % 逐代回溯：对每代 α 狼最优位置重新计算四指标，用于画各指标变化曲线
    for g = 1:Max_iter
        [~, mg] = mountain_objective(posHist(g,:), cfg);
        IndG(r,g,:) = [mg.Jr, mg.Jm, mg.Js, mg.Jl];
        KmatG(r,g)  = dot(cfg.eta, [mg.Jr, mg.Jm, mg.Js, mg.Jl]);
    end
    if s < gwoBest, gwoBest = s;  gwoX = pos; end
end
[~, gwoM] = mountain_objective(gwoX, cfg);
fprintf('GWO：Best=%.4f  Mean=%.4f  Worst=%.4f  Std=%.4e  平均耗时=%.1f ms\n', ...
    min(gwoScore), mean(gwoScore), max(gwoScore), std(gwoScore), mean(gwoTime));

%% ==================== 改进 IGWO（免疫克隆变异） ====================
igwoScore = zeros(runs,1);  igwoTime = zeros(runs,1);
igwoCurve = zeros(runs, Max_iter);
KmatI     = zeros(runs, Max_iter);
IndI      = zeros(runs, Max_iter, 4);
igwoBest = inf;  igwoX = [];  igwoM = [];  trigCnt = 0;

for r = 1:runs
    rng(r);                                 % 与 GWO 使用同一随机种子（配对实验）
    t0 = tic;
    [s, pos, conv, trig, posHist] = IGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, cfg);
    igwoTime(r) = toc(t0) * 1000;
    igwoScore(r) = s;  igwoCurve(r,:) = conv;  trigCnt = trigCnt + sum(trig);
    for g = 1:Max_iter
        [~, mg] = mountain_objective(posHist(g,:), cfg);
        IndI(r,g,:) = [mg.Jr, mg.Jm, mg.Js, mg.Jl];
        KmatI(r,g)  = dot(cfg.eta, [mg.Jr, mg.Jm, mg.Js, mg.Jl]);
    end
    if s < igwoBest, igwoBest = s;  igwoX = pos; end
end
[~, igwoM] = mountain_objective(igwoX, cfg);
fprintf('IGWO：Best=%.4f  Mean=%.4f  Worst=%.4f  Std=%.4e  平均耗时=%.1f ms\n', ...
    min(igwoScore), mean(igwoScore), max(igwoScore), std(igwoScore), mean(igwoTime));
fprintf('克隆变异累计触发次数（10 次运行合计）：%d\n\n', trigCnt);

%% ==================== 表1：无变异情况下的仿真结果 ====================
fprintf('===== 表1 无变异情况下的仿真结果 =====\n');
fprintf('%-12s', '序号');   for r = 0:runs-1, fprintf('%9d', r); end;  fprintf('%9s\n', '均值');
fprintf('%-12s', '时间/ms');  fprintf('%9.1f', gwoTime);                fprintf('%9.1f\n', mean(gwoTime));
fprintf('%-12s', '最终目标值'); fprintf('%9.4f', KmatG(:,end)');         fprintf('%9.4f\n', mean(KmatG(:,end)));
fprintf('%-12s', '距离指标');  fprintf('%9.4f', IndG(:,end,4)');         fprintf('%9.4f\n', mean(IndG(:,end,4)));
fprintf('%-12s', '安全指标');  fprintf('%9.4f', IndG(:,end,3)');         fprintf('%9.4f\n', mean(IndG(:,end,3)));
fprintf('%-12s', '隐身指标');  fprintf('%9.4f', IndG(:,end,1)');         fprintf('%9.4f\n', mean(IndG(:,end,1)));
fprintf('%-12s', '机动性指标'); fprintf('%9.4f', IndG(:,end,2)');        fprintf('%9.4f\n', mean(IndG(:,end,2)));

%% ==================== 表2：有变异情况下的仿真结果 ====================
fprintf('\n===== 表2 有变异情况下的仿真结果 =====\n');
fprintf('%-12s', '序号');   for r = 0:runs-1, fprintf('%9d', r); end;  fprintf('%9s\n', '均值');
fprintf('%-12s', '时间/ms');  fprintf('%9.1f', igwoTime);               fprintf('%9.1f\n', mean(igwoTime));
fprintf('%-12s', '最终目标值'); fprintf('%9.4f', KmatI(:,end)');         fprintf('%9.4f\n', mean(KmatI(:,end)));
fprintf('%-12s', '距离指标');  fprintf('%9.4f', IndI(:,end,4)');         fprintf('%9.4f\n', mean(IndI(:,end,4)));
fprintf('%-12s', '安全指标');  fprintf('%9.4f', IndI(:,end,3)');         fprintf('%9.4f\n', mean(IndI(:,end,3)));
fprintf('%-12s', '隐身指标');  fprintf('%9.4f', IndI(:,end,1)');         fprintf('%9.4f\n', mean(IndI(:,end,1)));
fprintf('%-12s', '机动性指标'); fprintf('%9.4f', IndI(:,end,2)');        fprintf('%9.4f\n', mean(IndI(:,end,2)));

%% ==================== 图4：总指标变化曲线（对应论文图4） ====================
figure('Name', '图4 总指标变化曲线', 'Color', 'w');
plot(1:Max_iter, mean(KmatI, 1), 'LineWidth', 1.6, 'Color', [0 0.45 0.74]);
xlabel('迭代次数'); ylabel('指标值'); title('总指标变化曲线'); grid on;

%% ==================== 图5：各指标变化曲线（对应论文图5） ====================
figure('Name', '图5 各指标变化曲线', 'Color', 'w');
plot(1:Max_iter, mean(IndI(:,:,4), 1), 'LineWidth', 1.4, 'Color', [0.2 0.55 0.2]); hold on;  % 距离
plot(1:Max_iter, mean(IndI(:,:,3), 1), 'LineWidth', 1.4, 'Color', [0.85 0.33 0.1]);          % 安全
plot(1:Max_iter, mean(IndI(:,:,1), 1), 'LineWidth', 1.4, 'Color', [0.49 0.18 0.56]);         % 隐身
plot(1:Max_iter, mean(IndI(:,:,2), 1), 'LineWidth', 1.4, 'Color', [0 0.45 0.74]);            % 机动
xlabel('迭代次数'); ylabel('指标值'); title('各指标变化曲线');
legend('距离', '安全', '隐身', '机动性', 'Location', 'best'); grid on;

%% ==================== 图6/图7：无变异 vs 有变异 总指标变化曲线 ====================
figure('Name', '图6/图7 总指标变化曲线对比', 'Color', 'w');
col = lines(runs);
subplot(1, 2, 1); hold on;
for r = 1:runs
    plot(1:Max_iter, KmatG(r,:), 'Color', col(r,:), 'LineWidth', 0.8);
end
xlabel('迭代次数'); ylabel('指标值'); title('图6 无变异条件下的总指标变化曲线'); grid on;
subplot(1, 2, 2); hold on;
for r = 1:runs
    plot(1:Max_iter, KmatI(r,:), 'Color', col(r,:), 'LineWidth', 0.8);
end
xlabel('迭代次数'); ylabel('指标值'); title('图7 有变异条件下的总指标变化曲线'); grid on;
lgd = arrayfun(@(k) sprintf('第%d次', k), 0:runs-1, 'UniformOutput', false);
legend(lgd, 'Location', 'best');

%% ==================== 结果对比 ====================
fprintf('\n===== 指标对比（各自最优航路） =====\n');
fprintf('%-10s %10s %10s %10s %10s %10s\n', '算法', '总指标K', '隐蔽Jr', '机动Jm', '安全Js', '航程Jl');
fprintf('%-10s %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'GWO', ...
    gwoBest, gwoM.Jr, gwoM.Jm, gwoM.Js, gwoM.Jl);
fprintf('%-10s %10.4f %10.4f %10.4f %10.4f %10.4f\n', 'IGWO', ...
    igwoBest, igwoM.Jr, igwoM.Jm, igwoM.Js, igwoM.Jl);

fprintf('\n航程 %.1f m（直线距离 %.1f m，航程余量 %.1f%%）\n', ...
    igwoM.Ltot, norm(cfg.T - cfg.S), 100*(igwoM.Ltot/norm(cfg.T-cfg.S) - 1));
fprintf('最小离地高度 %.1f m（安全余量要求 %.1f m）\n', min(igwoM.agl), cfg.h_safe);
fprintf('最大转角 %.1f 度（约束上限 %.1f 度）\n', max(igwoM.theta), cfg.theta_max);
if igwoM.feasible
    fprintf('约束校验：全部满足\n');
else
    fprintf('约束校验：存在违反，罚项 = %.4f\n', igwoM.pen);
end

fprintf('\n===== 克隆变异机制消融（论文表2/表3） =====\n');
fprintf('终值离散度(Std)：GWO=%.4e  →  IGWO=%.4e\n', std(gwoScore), std(igwoScore));
fprintf('平均耗时：GWO=%.1f ms  →  IGWO=%.1f ms  (增加 %.1f%%)\n', ...
    mean(gwoTime), mean(igwoTime), 100*(mean(igwoTime)/mean(gwoTime) - 1));
fprintf('综合指标提升：%.2f%%\n', 100*(mean(gwoScore) - mean(igwoScore))/mean(gwoScore));

%% ==================== 绘图：航路对比 ====================
cmp = struct('P', {gwoM.P, igwoM.P}, 'name', {'GWO（无变异）', 'IGWO（有变异）'});
plot_mountain(igwoM, cfg, '山地航路规划结果对比', cmp);

%% ==================== 绘图：收敛曲线对比 ====================
figure('Name', '收敛曲线对比', 'Color', 'w');
semilogy(mean(gwoCurve,1), 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]); hold on;
semilogy(mean(igwoCurve,1), 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xlabel('迭代次数'); ylabel('总指标 K（对数坐标）');
legend('GWO（无变异）', 'IGWO（有变异）', 'Location', 'northeast');
title(sprintf('收敛曲线对比（%d 次平均）', runs)); grid on;
