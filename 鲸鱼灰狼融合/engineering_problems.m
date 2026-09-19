function [lb, ub, dim, fobj, info] = engineering_problems(name)
% 返回三个工程问题的边界、维度、目标函数句柄与参考信息
% 输入：name 问题名('spring'/'welded'/'vessel')

PENALTY = 1e9;   % 罚系数，取值需远大于目标函数量级，确保可行解优先

switch lower(name)

    % ---------- 弹簧设计：3 个变量 ----------
    case 'spring'
        dim = 3;
        lb = [0.05, 0.25, 2];               % 变量下界：d, D, N
        ub = [2, 1.3, 15];                  % 变量上界：d, D, N
        info.name = '拉伸/压缩弹簧设计';
        info.x0 = [0.051690, 0.356750, 11.287126]; % 文献已知最优解
        info.f0 = 0.012665;                 % 文献已知最优值

    % ---------- 焊接梁设计：4 个变量 ----------
    case 'welded'
        dim = 4;
        lb = [0.1, 0.1, 0.1, 0.1];          % 变量下界：h, l, t, b
        ub = [2, 10, 10, 2];                % 变量上界：h, l, t, b
        info.name = '焊接梁设计';
        info.x0 = [0.205730, 3.470489, 9.036624, 0.205730];
        info.f0 = 1.724852;

    % ---------- 压力容器设计：4 个变量 ----------
    case 'vessel'
        dim = 4;
        lb = [0, 0, 10, 10];                % 变量下界：Ts, Th, R, L
        ub = [99, 99, 200, 200];            % 变量上界：Ts, Th, R, L
        info.name = '压力容器设计';
        info.x0 = [0.8125, 0.4375, 42.0984, 176.6366];
        info.f0 = 6059.7143;

    otherwise
        error('未知的工程问题：%s', name);
end

% 把罚系数固定进匿名函数，作为优化算法的目标函数
fobj = @(x) penalized_objective(name, x, PENALTY);
end