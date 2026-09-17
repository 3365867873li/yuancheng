function [K, m] = mountain_objective(x, cfg)
% ============================================================
% 计算山地航路的四项指标与加权总指标（论文 1.3 节，式(5)~(12)）
% 输入：x    决策变量 1×(3n)，依次为 [x1,y1,z1, x2,y2,z2, ...]
%       cfg  场景配置
% 输出：K    加权总指标（含约束罚项，越小越好）
%       m    结构体，保存四项指标与航路点，便于分析与绘图
% ============================================================

%% ---------- 第 1 步：解码航路 ----------
n = cfg.n;
W  = reshape(x, 3, n)';             % n×3，每行一个航路点
P  = [cfg.S; W; cfg.T];             % 拼接起点与终点
Np = size(P, 1);                    % 航路点总数 N_p

seg  = P(2:end, :) - P(1:end-1, :); % 各航段向量
Lseg = sqrt(sum(seg.^2, 2));        % 各航段长度 l_n
Ltot = sum(Lseg);                   % 总航程

%% ---------- 第 2 步：沿航路采样（用于地形与安全评估） ----------
ns   = cfg.samples_per_seg;
Pt   = zeros((Np-1)*ns, 3);         % 采样点坐标
sw   = zeros((Np-1)*ns, 1);         % 每个采样点代表的航路长度
k = 0;
for i = 1:Np-1
    for s = 1:ns
        k = k + 1;
        Pt(k, :) = P(i, :) + (s/ns) * seg(i, :);
        sw(k)    = Lseg(i) / ns;
    end
end
he_s = terrain_height(Pt(:,1), Pt(:,2), cfg);   % 采样点处的地面高度

%% ---------- 第 3 步：隐蔽性指标 J_r（式(6)(7)） ----------
he = terrain_height(P(:,1), P(:,2), cfg);       % 各航路点的地面高度
Jr_node = zeros(Np, 1);
for i = 1:Np
    if P(i,3) > cfg.h_ctrl
        Jr_node(i) = 0.2;                       % 超过控制高度：给定惩罚
    else
        % 离地高度越大越不隐蔽；下限截断为 0（低于地面由安全性指标处理）
        Jr_node(i) = max(0, P(i,3) - he(i)) / (200 * Np);
    end
end
Jr = cfg.k_r * sum(Jr_node);                    % K_r

%% ---------- 第 4 步：机动性指标 J_m（式(7)(8)） ----------
% 各节点处相邻两航段的夹角 θ_n（端点无转角，置 0）
theta = zeros(Np, 1);
for i = 2:Np-1
    v1 = seg(i-1, :);  v2 = seg(i, :);
    c  = dot(v1, v2) / (norm(v1) * norm(v2) + eps);
    theta(i) = acos(min(1, max(-1, c))) * 180 / pi;
end

th_set = 10;                                    % 设定转弯角 / 度
Jm_node = zeros(Np, 1);
for i = 1:Np
    th = theta(i);
    if th <= th_set
        Jm_node(i) = 0.4 * th / th_set;          % 小转角：线性惩罚
    elseif th <= cfg.theta_max
        Jm_node(i) = 0.4 + (th - th_set) / cfg.theta_max * (10 / Np);
    else                                        % 超过最大夹角：继续放大
        Jm_node(i) = 0.4 + (10 / Np) + (th - cfg.theta_max) / 10 * (10 / Np);
    end
end
Jm = cfg.k_m * sum(Jm_node);                    % K_m

%% ---------- 第 5 步：安全性指标 J_s（式(9)(10)） ----------
% h_min 为该点地面高度 + 安全余量，为动态值
h_min = he_s + cfg.h_safe;
Js_s  = zeros(size(Pt,1), 1);
below = Pt(:,3) < h_min;                        % 低于最低飞行高度的点
Js_s(below) = (h_min(below) - Pt(below,3)) / (100 * Np);
Js = sum(Js_s .* sw) / Ltot;                    % K_s（按航程加权平均）

%% ---------- 第 6 步：航程指标 J_l（式(11)） ----------
Jl = cfg.k_l * sum(Lseg);                       % K_l = k_l · Σ l_n

%% ---------- 第 7 步：约束罚项（论文 1.1.2 节的四类约束） ----------
pen = 0;
% (1) 最大航程约束：Σ l_n ≤ S_max               式(1)
pen = pen + max(0, Ltot - cfg.S_max) / cfg.S_max;
% (2) 航线夹角约束：θ ≤ θ_max                    式(3)
pen = pen + sum(max(0, (theta - cfg.theta_max) / cfg.theta_max));
% (3) 最长/最短航路点间隔约束：d_min ≤ l_n ≤ d_max 式(4)
pen = pen + sum(max(0, (Lseg - cfg.d_max) / cfg.d_max)) ...
          + sum(max(0, (cfg.d_min - Lseg) / cfg.d_min));
% (4) 飞行高度约束：h_max ≥ h_n ≥ h_min           式(5)
pen = pen + sum(max(0, (P(:,3) - cfg.h_max) / cfg.h_max));
% (5) 撞地约束：航路上任意点都不得低于地面
pen = pen + mean(max(0, he_s - Pt(:,3))) / 100;   % 平均触地深度归一化
% (6) 安全性指标本身作为软约束（保证进入可行域）
pen = pen + Js;

%% ---------- 第 8 步：加权求和得到总指标（式(12)） ----------
eta = cfg.eta;
K = eta(1)*Jr + eta(2)*Jm + eta(3)*Js + eta(4)*Jl + cfg.PENALTY * pen;

%% ---------- 第 9 步：打包输出 ----------
m.Jr = Jr;  m.Jm = Jm;  m.Js = Js;  m.Jl = Jl;
m.theta = theta;            % 各节点转角 / 度
m.Ltot  = Ltot;             % 总航程 / m
m.P     = P;                % 航路点序列
m.he    = he;               % 航路点处地面高度
m.agl   = P(:,3) - he;      % 离地高度
m.pen   = pen;              % 约束违反量
m.feasible = (pen <= 1e-9);

end

% ============================================================
function h = terrain_height(x, y, cfg)
% 由高斯峰体叠加得到地面高度 h(x,y)
h = zeros(numel(x), 1);
for i = 1:size(cfg.peaks, 1)
    px = cfg.peaks(i,1);  py = cfg.peaks(i,2);
    ph = cfg.peaks(i,3);  ps = cfg.peaks(i,4);
    h = h + ph * exp(-((x - px).^2 + (y - py).^2) / (2 * ps^2));
end
end
