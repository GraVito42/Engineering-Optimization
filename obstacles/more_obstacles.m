
% --- 1. Global Configuration ---
cfg.A = [0; 0]; cfg.B = [10; 0];

% MATRIX OF OBSTACLES: [Center_X, Center_Y, Radius]
% You can add as many rows as you want here
cfg.obstacles = [ 3,  1.0, 1.2; 
                  5, -1.5, 0.9;
                  7,  0.5, 0.9 
                  1,  0.4, 0.6
                  1,  -0.4, 0.2];
              
cfg.dt = 0.01; cfg.t = (0:cfg.dt:1)';

% Bernstein Basis (Degree 3)
cfg.b1 = 3 .* cfg.t .* (1 - cfg.t).^2;
cfg.b2 = 3 .* cfg.t.^2 .* (1 - cfg.t);
dir_vec = cfg.B - cfg.A;
cfg.n_vec = [-dir_vec(2); dir_vec(1)] / norm(dir_vec);

% --- 2. Grid Search with Multi-Collision Check ---
res = 1000; 
alpha_range = linspace(-15, 15, res);
[A1, A2] = meshgrid(alpha_range, alpha_range);

J_pure = zeros(res);
J_admissible = nan(res); 

for i = 1:res
    for j = 1:res
        % Generate curve
        P = (cfg.A' + cfg.t*(cfg.B - cfg.A)') + (A1(i,j)*cfg.b1 + A2(i,j)*cfg.b2)*cfg.n_vec';
        
        % Compute Base Cost
        cost = bernstein_cost(P, cfg.dt);
        J_pure(i,j) = cost;
        
        % Check ALL obstacles
        is_safe = true;
        for k = 1:size(cfg.obstacles, 1)
            obs_C = cfg.obstacles(k, 1:2);
            obs_R = cfg.obstacles(k, 3);
            if obstacle_distance(obs_C, obs_R, P) < 0
                is_safe = false;
                break; % Exit loop if one collision is found
            end
        end
        
        if is_safe
            J_admissible(i,j) = cost;
        end
    end
end

% Find numerical optimum
[min_cost, idx] = min(J_admissible(:), [], 'omitnan');
if isempty(idx) || isnan(min_cost)
    error('No admissible path found. Try increasing alpha range or reducing obstacles.');
end
[r, c] = ind2sub(size(J_admissible), idx);
opt_alpha = [A1(r,c), A2(r,c)];

% --- 3. Visualization ---

% PLOT 1: 3D Surface with Multiple Holes
figure('Name', 'Multi-Obstacle Energy Landscape', 'Color', 'w'); hold on;
surf(A1, A2, J_pure, 'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
surf(A1, A2, J_admissible, 'FaceColor', 'r', 'FaceAlpha', 0.8, 'EdgeColor', 'none');
plot3(opt_alpha(1), opt_alpha(2), min_cost, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
title('3D Surface with Multiple Constraints (Red = Admissible)');
view(-45, 35); camlight; lighting phong; grid on;

% PLOT 2: 2D Multi-Obstacle Path
figure('Name', 'Multi-Obstacle Optimal Path', 'Color', 'w'); hold on; axis equal; grid on;

% Draw All Obstacles
th = linspace(0, 2*pi, 100);
for k = 1:size(cfg.obstacles, 1)
    ox = cfg.obstacles(k,1); oy = cfg.obstacles(k,2); or = cfg.obstacles(k,3);
    fill(ox + or*cos(th), oy + or*sin(th), 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'r');
end

% Re-generate and Plot Optimal Path
P_opt = (cfg.A' + cfg.t*(cfg.B - cfg.A)') + (opt_alpha(1)*cfg.b1 + opt_alpha(2)*cfg.b2)*cfg.n_vec';
plot(P_opt(:,1), P_opt(:,2), 'b', 'LineWidth', 3);
plot([cfg.A(1), cfg.B(1)], [cfg.A(2), cfg.B(2)], 'ko', 'MarkerFaceColor', 'k');
title(sprintf('Optimal Path: \\alpha_1=%.1f, \\alpha_2=%.1f', opt_alpha(1), opt_alpha(2)));
