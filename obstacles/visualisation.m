clc
clear all
close all

% --- Global Configuration ---
cfg.A = [0; 0]; cfg.B = [10; 2];
cfg.obs_C = [5; 1.5]; cfg.obs_R = 1.2;
cfg.dt = 0.01; cfg.t = (0:cfg.dt:1)';

% Bernstein Basis (Degree 3)
cfg.b1 = 3 .* cfg.t .* (1 - cfg.t).^2;
cfg.b2 = 3 .* cfg.t.^2 .* (1 - cfg.t);
dir_vec = cfg.B - cfg.A;
cfg.n_vec = [-dir_vec(2); dir_vec(1)] / norm(dir_vec);

% --- 1. Grid Search for Optimum ---
res = 60; % Higher resolution for better optimum accuracy
alpha_range = linspace(-15, 15, res);
[A1, A2] = meshgrid(alpha_range, alpha_range);

J_pure = zeros(res);
J_admissible = nan(res); % Initialize with NaN

for i = 1:res
    for j = 1:res
        % Generate curve for this specific alpha pair
        P = (cfg.A' + cfg.t*(cfg.B - cfg.A)') + (A1(i,j)*cfg.b1 + A2(i,j)*cfg.b2)*cfg.n_vec';
        
        % Compute Cost and Distance
        cost = bernstein_cost(P, cfg.dt);
        J_pure(i,j) = cost;
        
        d = obstacle_distance(cfg.obs_C, cfg.obs_R, P);
        if d >= 0
            J_admissible(i,j) = cost;
        end
    end
end

% Find the numerical optimum in the admissible set
[min_cost, idx] = min(J_admissible(:), [], 'omitnan');
[r, c] = ind2sub(size(J_admissible), idx);
opt_alpha = [A1(r,c), A2(r,c)];

% --- 2. Visualization ---

% PLOT 1: 3D Surface
figure('Name', '3D Energy Landscape', 'Color', 'w'); hold on;
surf(A1, A2, J_pure, 'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
surf(A1, A2, J_admissible, 'FaceColor', 'r', 'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot3(opt_alpha(1), opt_alpha(2), min_cost, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
title(sprintf('Optimum found at \\alpha_1=%.2f, \\alpha_2=%.2f', opt_alpha(1), opt_alpha(2)));
view(-45, 35); camlight; lighting phong; grid on;

% PLOT 2: 2D Optimal Path
figure('Name', 'Optimal Trajectory Result', 'Color', 'w'); hold on; axis equal; grid on;

% Draw Obstacle
th = linspace(0, 2*pi, 100);
fill(cfg.obs_C(1) + cfg.obs_R*cos(th), cfg.obs_C(2) + cfg.obs_R*sin(th), ...
     'r', 'FaceAlpha', 0.3, 'EdgeColor', 'r');

% Re-generate Optimal Path
P_opt = (cfg.A' + cfg.t*(cfg.B - cfg.A)') + (opt_alpha(1)*cfg.b1 + opt_alpha(2)*cfg.b2)*cfg.n_vec';

% Plot start/end and path
plot(P_opt(:,1), P_opt(:,2), 'b', 'LineWidth', 3);
plot([cfg.A(1), cfg.B(1)], [cfg.A(2), cfg.B(2)], 'ko', 'MarkerFaceColor', 'k');

title('Optimal Path (Minimum Energy + Avoidance)');
xlabel('X [m]'); ylabel('Y [m]');
legend('Obstacle', 'Optimal Path (Admissible)');