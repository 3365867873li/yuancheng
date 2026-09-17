clear; clc; close all;   % 清空变量、命令行、关闭旧图窗

%% ==================== 用户控制区：选择要运行的工程问题 ====================
% 可选值：'spring'(弹簧) / 'welded'(焊接梁) / 'vessel'(压力容器) / 'all'(全部)
% 也可写成元胞数组同时选多个，例如：{'spring','vessel'}
PROBLEM_SELECT = 'spring';   % <<< 只需修改这一行即可控制运行哪个问题

USE_MENU = true;            % <<< 设为 true 则弹出菜单让你交互勾选（支持多选）
                             %     设为 true 时，上面的 PROBLEM_SELECT 被忽略

%% ==================== 参数设置(与论文一致) ====================
SearchAgents_no = 30;    % 狼群规模
Max_iter = 500;          % 最大迭代次数
runs = 30;               % 独立重复运行次数，用于统计鲁棒性

% 全部可用问题的名称与中文标签
allNames  = {'spring', 'welded', 'vessel'};
allLabels = {'拉伸/压缩弹簧设计', '焊接梁设计', '压力容器设计'};

%% ==================== 解析用户选择：确定本次要运行哪些问题 ====================
if USE_MENU
    % --- 方式二：弹出列表菜单，交互多选 ---
    [idx, ok] = listdlg('PromptString', '选择要运行的工程问题（可按住 Ctrl 多选）', ...
        'SelectionMode', 'multiple', ...
        'ListString', allLabels, ...
        'ListSize', [260 120], ...
        'InitialValue', 1:numel(allLabels));
    if ~ok
        disp('已取消选择，程序退出。');   % 用户点关闭或取消
        return;
    end
    problems = allNames(idx);            % 取出对应的英文问题名
else
    % --- 方式一：使用脚本顶部变量控制 ---
    if ischar(PROBLEM_SELECT) && strcmpi(PROBLEM_SELECT, 'all')
        problems = allNames;             % 全部问题
    elseif ischar(PROBLEM_SELECT)
        problems = {PROBLEM_SELECT};     % 单个问题
    else
        problems = PROBLEM_SELECT;       % 元胞数组，多个问题
    end

    % 统一转小写，便于大小写不敏感匹配
    problems = cellfun(@lower, problems, 'UniformOutput', false);

    % 校验名称合法性，防止手误
    for k = 1:numel(problems)
        if ~ismember(problems{k}, allNames)
            error('无效的问题名："%s"（可选 spring / welded / vessel / all）', problems{k});
        end
    end
end

fprintf('本次将运行 %d 个工程问题：%s\n', numel(problems), strjoin(problems, ', '));

%% ==================== 逐个问题求解 ====================
for p = 1:numel(problems)
    name = problems{p};
    % 读取该问题的上下界、维度、目标函数与参考信息
    [lb, ub, dim, fobj, info] = engineering_problems(name);

    scores    = zeros(runs, 1);          % 存放每次运行的最优适应度
    positions = zeros(runs, dim);        % 存放每次运行的最优解
    curves    = zeros(runs, Max_iter);   % 存放每次运行的收敛曲线

    % 独立运行 runs 次，消除随机性影响
    for r = 1:runs
        [s, pos, conv] = GWO(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
        scores(r)    = s;                % 记录最优适应度
        positions(r, :) = pos;           % 记录最优解
        curves(r, :) = conv;             % 记录收敛曲线
    end

    % 找出所有运行中最优的一次
    [~, idx] = min(scores);
    bestX = positions(idx, :);
    % 用真实目标(不含罚项)与约束评估该最优解，便于与文献对比
    [fTrue, gBest] = engineering_eval(name, bestX);

    % ---------------- 打印结果 ----------------
    fprintf('\n===== %s =====\n', info.name);
    fprintf('最优设计变量 x*   : %s\n', mat2str(bestX, 6));
    fprintf('最优目标值        : %.6f\n', fTrue);
    fprintf('文献参考最优值    : %.6f   x_ref = %s\n', info.f0, mat2str(info.x0, 6));
    fprintf('最大约束违反量    : %.3e  (<=0 表示可行)\n', max(gBest));
    fprintf('--- %d 次独立运行统计 ---\n', runs);
    fprintf('Best  = %.6f\n', min(scores));   % 最优
    fprintf('Mean  = %.6f\n', mean(scores));  % 平均
    fprintf('Worst = %.6f\n', max(scores));   % 最差
    fprintf('Std   = %.6e\n', std(scores));   % 标准差(鲁棒性)

    % ---------------- 画平均收敛曲线 ----------------
    figure('Name', info.name);
    semilogy(mean(curves, 1), 'LineWidth', 1.5);  % 对数纵轴显示平均收敛
    xlabel('迭代次数'); ylabel('最优适应度(对数坐标)');
    title([info.name ' - 30 次运行平均收敛曲线']);
    grid on;                                      % 打开网格
end