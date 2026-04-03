clc
clear all
close all

% --- 1. SETTINGS & CONFIGURATION ---
cfg.n_vars = 5; % Number of variables
cfg.A = [0; 0]; cfg.B = [10; 0];
cfg.obstacles = [ 3,  1.0, 1.2; 
                  5, -1.5, 0.9;
                  7,  0.5, 0.9 
                  1,  0.4, 0.3
                  1,  -0.4, 0.2];
cfg.dt = 0.01; cfg.t = (0:cfg.dt:1)';

dir_vec = cfg.B - cfg.A;
cfg.n_vec = [-dir_vec(2); dir_vec(1)] / norm(dir_vec);

% --- 2. DYNAMIC BASIS GENERATION ---
d = cfg.n_vars + 1; 
cfg.M = zeros(length(cfg.t), cfg.n_vars);
for i = 1:cfg.n_vars
    cfg.M(:, i) = nchoosek(d, i) .* (cfg.t.^i) .* ((1-cfg.t).^(d-i));
end

% --- 3. OPTIMIZATION ---
opt_alphas = zeros(1, cfg.n_vars);
if cfg.n_vars == 2
    res = 1000; vals = linspace(-12, 12, res);
    [A1, A2] = meshgrid(vals, vals);
    J_pure = zeros(res); J_adm = nan(res);
    for i = 1:res
        for j = 1:res
            [cost, safe] = evaluate_config([A1(i,j), A2(i,j)], cfg);
            J_pure(i,j) = cost;
            if safe, J_adm(i,j) = cost; end
        end
    end
    [min_cost, idx] = min(J_adm(:), [], 'omitnan');
    if ~isempty(idx)
        [r, c] = ind2sub(size(J_adm), idx);
        opt_alphas = [A1(r,c), A2(r,c)];
        plot_3d_surface(A1, A2, J_pure, J_adm, opt_alphas, min_cost);
    end
end

launch_interactive_2d(cfg, opt_alphas);


% --- UPDATED INTERACTIVE PLOTTER (Fixes the Dot Indexing Error) ---
function launch_interactive_2d(cfg, initial_alphas)
    fig = figure('Name', 'Interactive Planner', 'Color', 'w', 'Position', [100 100 900 650]);
    ax = axes('Position', [0.1 0.4 0.8 0.5]); hold on; grid on; axis equal;
    
    % Draw Obstacles
    th = linspace(0, 2*pi, 50);
    for k = 1:size(cfg.obstacles, 1)
        fill(cfg.obstacles(k,1)+cfg.obstacles(k,3)*cos(th), ...
             cfg.obstacles(k,2)+cfg.obstacles(k,3)*sin(th), 'r', 'FaceAlpha', 0.2);
    end
    
    path_h = plot(0,0, 'b', 'LineWidth', 2.5);
    plot([cfg.A(1), cfg.B(1)], [cfg.A(2), cfg.B(2)], 'ko', 'MarkerFaceColor', 'k');
    
    % Create N Sliders using old-school handles for compatibility
    h_sliders = zeros(1, cfg.n_vars);
    for i = 1:cfg.n_vars
        y_pos = 20 + (i-1)*(250/cfg.n_vars);
        uicontrol('Style', 'text', 'Position', [50 y_pos 50 20], 'String', sprintf('A%d', i));
        
        % Store handle in array
        h_sliders(i) = uicontrol('Style', 'slider', 'Min', -15, 'Max', 15, ...
                                 'Value', initial_alphas(i), ...
                                 'Position', [110 y_pos 600 20]);
    end

    % Set the callback for all sliders at once
    % Using 'get' inside the callback avoids the dot indexing error
    for i = 1:cfg.n_vars
        set(h_sliders(i), 'Callback', @(src, ev) update_ui());
    end

    update_ui();

    function update_ui()
        % CORRECT WAY to access values from numeric handles
        alphas = zeros(1, cfg.n_vars);
        for k = 1:cfg.n_vars
            alphas(k) = get(h_sliders(k), 'Value');
        end
        
        [cost, safe] = evaluate_config(alphas, cfg);
        P = (cfg.A' + cfg.t*(cfg.B - cfg.A)') + (cfg.M * alphas') * cfg.n_vec';
        
        color = 'b'; if ~safe, color = 'r'; end
        set(path_h, 'XData', P(:,1), 'YData', P(:,2), 'Color', color);
        title(ax, sprintf('Total Cost: %.2f | Safe: %d', cost, safe));
    end
end

% --- SUPPORT FUNCTIONS ---
function [cost, is_safe] = evaluate_config(alphas, cfg)
    P = (cfg.A' + cfg.t*(cfg.B - cfg.A)') + (cfg.M * alphas') * cfg.n_vec';
    v = diff(P)/cfg.dt; a = diff(v)/cfg.dt;
    cost = sum(sqrt(sum(v.^2, 2)))*cfg.dt + 0.5*sum(sum(a.^2, 2))*cfg.dt;
    is_safe = true;
    for k = 1:size(cfg.obstacles, 1)
        dist = min(sqrt((P(:,1)-cfg.obstacles(k,1)).^2 + (P(:,2)-cfg.obstacles(k,2)).^2)) - cfg.obstacles(k,3);
        if dist < 0, is_safe = false; break; end
    end
end

function plot_3d_surface(A1, A2, Jp, Ja, opt, m_val)
    figure('Name', 'Cost Landscape', 'Color', 'w'); hold on;
    surf(A1, A2, Jp, 'FaceColor', 'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    surf(A1, A2, Ja, 'FaceColor', 'r', 'FaceAlpha', 0.85, 'EdgeColor', 'none');
    plot3(opt(1), opt(2), m_val, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    view(-45, 35); camlight; lighting phong; grid on;
end