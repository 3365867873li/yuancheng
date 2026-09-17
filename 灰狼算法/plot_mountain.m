function plot_mountain(m, cfg, figTitle, cmp)
% ============================================================
% 绘制山地航路规划结果
% 输入：m         mountain_objective 返回的指标结构体（含 m.P）
%       cfg       场景配置
%       figTitle  图窗标题（可选）
%       cmp       可选：结构体数组 {P1,P2,...} 用于多算法航路对比
% ============================================================

if nargin < 3 || isempty(figTitle), figTitle = '山地航路规划'; end

%% ---------- 生成地形网格 ----------
gx = linspace(min(cfg.S(1), cfg.T(1)) - 800, max(cfg.S(1), cfg.T(1)) + 800, 80);
gy = linspace(cfg.y_range(1), cfg.y_range(2), 80);
[GX, GY] = meshgrid(gx, gy);
GZ = zeros(size(GX));
for i = 1:size(cfg.peaks, 1)
    px = cfg.peaks(i,1); py = cfg.peaks(i,2);
    ph = cfg.peaks(i,3); ps = cfg.peaks(i,4);
    GZ = GZ + ph * exp(-((GX-px).^2 + (GY-py).^2) / (2*ps^2));
end

figure('Name', figTitle, 'Color', 'w');

%% ---------- 子图 1：三维山地航路 ----------
subplot(1, 2, 1);
hold on; grid on; box on;

surf(GX, GY, GZ, 'EdgeColor', 'none', 'FaceAlpha', 0.85);   % 山地曲面
colormap(subplot(1,2,1), parula);
shading interp;

% 参考直线航路
plot3([cfg.S(1), cfg.T(1)], [cfg.S(2), cfg.T(2)], [cfg.S(3), cfg.T(3)], ...
    'k:', 'LineWidth', 1.2);

if nargin >= 4 && ~isempty(cmp)
    cols = lines(numel(cmp));
    for i = 1:numel(cmp)
        Q = cmp(i).P;
        plot3(Q(:,1), Q(:,2), Q(:,3), '-o', 'LineWidth', 2, ...
            'Color', cols(i,:), 'MarkerSize', 4, 'MarkerFaceColor', 'w');
    end
    legend([{'地形', '直线'}, {cmp.name}], 'Location', 'best');
else
    P = m.P;
    plot3(P(:,1), P(:,2), P(:,3), '-o', 'LineWidth', 2, ...
        'Color', [0, 0.35, 0.85], 'MarkerSize', 5, 'MarkerFaceColor', 'w');
end

plot3(cfg.S(1), cfg.S(2), cfg.S(3), 'p', 'MarkerSize', 15, ...
    'MarkerFaceColor', [0, 0.7, 0], 'MarkerEdgeColor', 'k');
plot3(cfg.T(1), cfg.T(2), cfg.T(3), '^', 'MarkerSize', 11, ...
    'MarkerFaceColor', [0.85, 0.1, 0.1], 'MarkerEdgeColor', 'k');

xlabel('x / m'); ylabel('y / m'); zlabel('z / m');
title('三维山地航路'); view(-40, 28);

%% ---------- 子图 2：航路离地高度（检验是否撞山） ----------
subplot(1, 2, 2);
hold on; grid on; box on;

P = m.P;
% 沿航路加密采样，画离地高度曲线
tt  = linspace(0, 1, 400);
PP  = interp1(1:size(P,1), P, 1 + tt*(size(P,1)-1), 'linear');
hh  = zeros(size(PP,1),1);
for i = 1:size(cfg.peaks,1)
    px = cfg.peaks(i,1); py = cfg.peaks(i,2);
    ph = cfg.peaks(i,3); ps = cfg.peaks(i,4);
    hh = hh + ph * exp(-((PP(:,1)-px).^2 + (PP(:,2)-py).^2) / (2*ps^2));
end
s = [0; cumsum(sqrt(sum(diff(PP).^2, 2)))];      % 沿航路的弧长
plot(s, hh, '-', 'LineWidth', 1.5, 'Color', [0.5 0.35 0.1]);   % 地形高度
plot(s, PP(:,3), '-', 'LineWidth', 2, 'Color', [0, 0.35, 0.85]); % 航路高度
plot(s, hh + cfg.h_safe, '--', 'LineWidth', 1, 'Color', [0.85 0.1 0.1]);
legend({'地形高度', '航路高度', '最低飞行高度'}, 'Location', 'best');
xlabel('沿航路距离 / m'); ylabel('高度 / m');
title('航路高度与地形剖面对比');

end
