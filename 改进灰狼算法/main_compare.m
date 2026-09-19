clear; clc; close all;

%% ==================== 路径与参数设置 ====================
% 本文件夹已自包含：GWO.m / IGWO.m / benchmark_functions.m / initialization.m
% 与上一实验的"灰狼算法"文件夹相互独立，互不依赖

SearchAgents_no = 30;     % 狼群规模
Max_iter        = 500;    % 最大迭代次数
runs            = 30;     % 独立运行次数（GWO/IGWO 使用同一随机种子，配对公平对比）
Functions       = 1:13;   % 测试 F1~F13

% 免疫克隆变异参数（与山地航路实验一致，论文式(16)，见 mountain_config.m）
cfg.immu.Nc     = 6;       % 克隆规模
cfg.immu.pm     = 0.25;    % 逐维变异概率
cfg.immu.sigma  = 0.08;    % 变异步长系数（相对变量范围）
cfg.immu.gap    = 4;       % 触发判断的代数间隔（式(16)：g-4）
cfg.immu.minGen = 30;      % 触发的最小代数（式(16)：g > 30）

% Lévy 飞行参数（IGWO 新增全局探索算子，可调）
cfg.levy.p      = 0.3;     % 每只狼每代施加 Lévy 扰动的概率
cfg.levy.lambda = 1.5;     % 重尾指数 λ∈(1,2]，越大跃迁越"重尾"

nF = numel(Functions);

%% ==================== 结果存储 ====================
gwoRes = struct('Best', zeros(nF,1), 'Mean', zeros(nF,1), ...
                'Worst', zeros(nF,1), 'Std', zeros(nF,1), 'Curve', {cell(nF,1)});
igwoRes = gwoRes;               % IGWO 的结果结构（同字段）
win     = repmat('=', nF, 1);   % 每函数胜者：I=IGWO  G=GWO  ==持平（按均值）
trigRate = zeros(nF, 1);        % IGWO 变异触发率（触发代数 / 总迭代）

%% ==================== 主循环：逐函数对比 ====================
for i = 1:nF
    F = Functions(i);
    [lb, ub, dim, fobj] = benchmark_functions(F);       % 取函数配置

    % ---- 触发阈值采用论文式(16)原值 eps=0.01（与山地航路实验一致） ----
    cfg.immu.eps = 0.01;

    % ---- 配对运行 GWO 与 IGWO ----
    gwoS = zeros(runs,1);  gwoC = zeros(runs, Max_iter);
    igS  = zeros(runs,1);  igC  = zeros(runs, Max_iter);
    trigCnt = zeros(runs,1);

    for r = 1:runs
        rng(r);                                             % 同一随机流，公平对比
        [s1, ~, c1] = GWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
        gwoS(r) = s1;  gwoC(r,:) = c1;

        rng(r);                                             % 复位到相同种子
        [s2, ~, c2, trig] = IGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj, cfg);
        igS(r) = s2;  igC(r,:) = c2;  trigCnt(r) = sum(trig);
    end

    % ---- 统计四项指标 ----
    gwoRes.Best(i)  = min(gwoS);   gwoRes.Worst(i) = max(gwoS);
    gwoRes.Mean(i)  = mean(gwoS);  gwoRes.Std(i)   = std(gwoS);
    gwoRes.Curve{i} = mean(gwoC, 1);

    igwoRes.Best(i)  = min(igS);   igwoRes.Worst(i) = max(igS);
    igwoRes.Mean(i)  = mean(igS);  igwoRes.Std(i)   = std(igS);
    igwoRes.Curve{i} = mean(igC, 1);

    trigRate(i) = mean(trigCnt) / Max_iter;

    % ---- 判断该函数上谁更优（均值 + 相对容差，避免把数值噪声当胜负） ----
    mG = gwoRes.Mean(i);  mI = igwoRes.Mean(i);
    tol = 1e-6 * max(1, min(abs(mG), abs(mI)));
    if mI < mG - tol
        win(i) = 'I';
    elseif mG < mI - tol
        win(i) = 'G';
    end
end

%% ==================== 打印结果表 ====================
fprintf('\n==================== 基准函数对比：GWO vs IGWO ====================\n');
fprintf('%-5s %-24s %-24s %-8s %-8s\n', '函数', 'GWO Mean', 'IGWO Mean', '胜者', '触发率');
for i = 1:nF
    fprintf('F%-4d %-24s %-24s %-8s %6.1f%%\n', Functions(i), ...
        sprintf('%.4e', gwoRes.Mean(i)), sprintf('%.4e', igwoRes.Mean(i)), win(i), 100*trigRate(i));
end

% ---- 汇总胜负计数 ----
nI = sum(win == 'I');  nG = sum(win == 'G');  nT = sum(win == '=');
fprintf('\n========== 汇总 ==========\n');
fprintf('IGWO 均值更优的函数：%d 个（', nI);
fprintf('%s', strjoin(arrayfun(@(i) sprintf('F%d', Functions(i)), find(win=='I'), 'UniformOutput', false), ', '));
fprintf('）\nGWO  均值更优的函数：%d 个（', nG);
fprintf('%s', strjoin(arrayfun(@(i) sprintf('F%d', Functions(i)), find(win=='G'), 'UniformOutput', false), ', '));
fprintf('）\n持平：%d 个\n', nT);

% 在 IGWO 胜出的函数上统计平均提升幅度
if nI > 0
    idx = find(win == 'I');
    imp = (gwoRes.Mean(idx) - igwoRes.Mean(idx)) ./ max(abs(gwoRes.Mean(idx)), 1e-300);
    fprintf('IGWO 胜出函数的平均提升幅度：%.2f%%\n', 100 * mean(imp));
end
if nG > 0
    idx = find(win == 'G');
    imp = (igwoRes.Mean(idx) - gwoRes.Mean(idx)) ./ max(abs(igwoRes.Mean(idx)), 1e-300);
    fprintf('GWO  胜出函数的平均领先幅度：%.2f%%\n', 100 * mean(imp));
end

%% ==================== 图1：各函数收敛曲线对比（4×4 子图） ====================
figure('Name', 'GWO vs IGWO 收敛曲线对比', 'Color', 'w');
for i = 1:nF
    subplot(4, 4, i);  hold on;
    cG = gwoRes.Curve{i};  cI = igwoRes.Curve{i};
    if min([cG cI]) > 0                            % 全为正：对数坐标更清晰
        semilogy(1:Max_iter, cG, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
        semilogy(1:Max_iter, cI, 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
    else                                          % 含负值（如 F8）：改用线性坐标
        plot(1:Max_iter, cG, 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
        plot(1:Max_iter, cI, 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
    end
    title(sprintf('F%d（胜者：%s）', Functions(i), win(i)));
    set(gca, 'FontSize', 8);  grid on;
end
subplot(4, 4, 14);  axis off;
subplot(4, 4, 15);  axis off;
subplot(4, 4, 16);  axis off;
subplot(4, 4, 13);
legend({'GWO', 'IGWO'}, 'Location', 'best', 'FontSize', 8);

%% ==================== 图2：均值对比柱状图（对数坐标） ====================
figure('Name', '均值对比柱状图', 'Color', 'w');
bar([max(abs(gwoRes.Mean), 1e-300), max(abs(igwoRes.Mean), 1e-300)], 'grouped');
set(gca, 'YScale', 'log');
xticklabels(arrayfun(@(i) sprintf('F%d', i), 1:nF, 'UniformOutput', false));
xlabel('基准函数'); ylabel('|Mean|（对数坐标）');
legend({'GWO', 'IGWO'}, 'Location', 'best'); grid on;
title('GWO vs IGWO 各函数最优值均值对比');
