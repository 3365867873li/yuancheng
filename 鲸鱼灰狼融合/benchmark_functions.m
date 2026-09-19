function [lb, ub, dim, fobj] = benchmark_functions(F)
% 返回第 F 个基准测试函数的边界、维度与函数句柄
% F1-F7 为单峰函数(测开发能力)，F8-F13 为多峰函数(测探索能力)

switch F
    case {1,2,3,4,5,6}      % 常见单峰函数搜索范围
        lb = -100; ub = 100; dim = 30;
    case 7                  % 噪声函数，范围小
        lb = -1.28; ub = 1.28; dim = 30;
    case 8                  % 多峰，范围大
        lb = -500; ub = 500; dim = 30;
    case 9                  % Rastrigin
        lb = -5.12; ub = 5.12; dim = 30;
    case 10                 % Ackley
        lb = -32; ub = 32; dim = 30;
    case 11                 % Griewank
        lb = -600; ub = 600; dim = 30;
    case {12,13}            % Penalized 系列
        lb = -50; ub = 50; dim = 30;
    otherwise
        error('F 必须为 1..13');
end
% 用匿名函数把函数编号与求值函数绑定
fobj = @(x) eval_F(F, x);
end

function o = eval_F(F, x)
% 根据编号计算具体测试函数值
dim = size(x, 2);           % 读取实际维度
switch F
    case 1
        o = sum(x.^2);                                          % 平方和
    case 2
        o = sum(abs(x)) + prod(abs(x));                         % 绝对值之和+乘积
    case 3
        o = sum(cumsum(x).^2);                                  % 累积平方和
    case 4
        o = max(abs(x));                                        % 最大值
    case 5
        o = sum(100*(x(2:dim) - x(1:dim-1).^2).^2 + (x(1:dim-1) - 1).^2); % Rosenbrock
    case 6
        o = sum((x + 0.5).^2);                                  % 取整平方和
    case 7
        o = sum((1:dim).*(x.^4)) + rand;                        % 四次方+噪声
    case 8
        o = sum(-x .* sin(sqrt(abs(x))));                       % Schwefel
    case 9
        o = sum(x.^2 - 10*cos(2*pi*x) + 10);                    % Rastrigin
    case 10
        o = -20*exp(-0.2*sqrt(sum(x.^2)/dim)) - exp(sum(cos(2*pi*x))/dim) + 20 + exp(1); % Ackley
    case 11
        o = sum(x.^2)/4000 - prod(cos(x./sqrt(1:dim))) + 1;     % Griewank
    case 12
        % Penalized 1：带正弦项与边界惩罚
        o = (pi/dim) * ( 10*(sin(pi*(1+(x(1)+1)/4)))^2 ...
            + sum((((x(1:dim-1)+1)/4).^2) .* (1 + 10*(sin(pi*(1+(x(2:dim)+1)/4))).^2)) ...
            + ((x(dim)+1)/4)^2 ) + sum(Ufun(x, 10, 100, 4));
    case 13
        % Penalized 2：另一种正弦惩罚形式
        o = 0.1 * ( (sin(3*pi*x(1)))^2 ...
            + sum((x(1:dim-1)-1).^2 .* (1 + (sin(3*pi*x(2:dim))).^2)) ...
            + (x(dim)-1)^2 * (1 + (sin(2*pi*x(dim)))^2) ) + sum(Ufun(x, 5, 100, 4));
end
end

function o = Ufun(x, a, k, m)
% 边界惩罚辅助函数：超出 [-a, a] 的部分按幂次 k 施加惩罚
o = k .* ((x - a).^m) .* (x > a) + k .* ((-x - a).^m) .* (x < -a);
end