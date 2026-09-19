clear; clc; close all;   % 清理环境

% ---------------- 参数设置（与论文一致） ----------------
SearchAgents_no = 30;    % 种群规模（论文：鲸鱼/狼群混合）
Max_iter = 500;          % 最大迭代次数
runs = 30;               % 独立运行次数（用于统计鲁棒性）
Functions = 1:13;        % 测试 F1 到 F13（F1-F7 单峰测开发，F8-F13 多峰测探索）

% 预分配统计结果数组
best_all  = zeros(1, numel(Functions));
mean_all  = zeros(1, numel(Functions));
worst_all = zeros(1, numel(Functions));
std_all   = zeros(1, numel(Functions));

% ---------------- 逐个函数求解并统计 ----------------
for i = 1:numel(Functions)
    F = Functions(i);
    [lb, ub, dim, fobj] = benchmark_functions(F);  % 取函数配置

    scores = zeros(runs, 1);          % 每次运行的最优值
    curves = zeros(runs, Max_iter);   % 每次运行的收敛曲线

    for r = 1:runs
        [s, ~, conv] = WOAGWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
        scores(r) = s;                % 记录本次最优值
        curves(r, :) = conv;          % 记录收敛曲线
    end

    % 统计四项指标
    best_all(i)  = min(scores);       % 最优
    worst_all(i) = max(scores);       % 最差
    mean_all(i)  = mean(scores);      % 平均
    std_all(i)   = std(scores);       % 标准差

    % 打印单行统计结果
    fprintf('F%-2d  Best=%.4e  Mean=%.4e  Worst=%.4e  Std=%.4e\n', ...
        F, best_all(i), mean_all(i), worst_all(i), std_all(i));

    % 画该函数的平均收敛曲线
    figure('Name', ['F' num2str(F)]);
    semilogy(mean(curves, 1), 'LineWidth', 1.5);
    xlabel('迭代次数'); ylabel('平均最优适应度');
    title(['基准函数 F' num2str(F) ' - 30 次运行平均收敛曲线']);
    grid on;
end

% ---------------- 汇总成表格 ----------------
disp(table(Functions', best_all', mean_all', worst_all', std_all', ...
    'VariableNames', {'F', 'Best', 'Mean', 'Worst', 'Std'}));
