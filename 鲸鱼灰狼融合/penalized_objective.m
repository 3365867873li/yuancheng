function f = penalized_objective(name, x, penalty)
% 用静态罚函数把约束问题转成无约束问题
% 输入：name 问题名；x 决策变量；penalty 罚系数
% 输出：f 惩罚后的适应度 = 原目标 + 罚项

% 先算出原目标 f0 和约束 g
[f0, g] = engineering_eval(name, x);

% 只对违反的约束(>0 的部分)施加惩罚，满足的约束不罚
% max(0, g) 把所有违反量提取出来，罚系数足够大即可逼迫解进入可行域
f = f0 + penalty * sum(max(0, g));
end