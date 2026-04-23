function stop = TrajectoryPlotter(x, optimValues, state, hfig)
    stop = false;
    global Par
    global details
    global ub lb
    global var_history

    n = length(x);

    zeroc = {};

    if strcmp(state, 'init') || isempty(var_history)
        var_history = table(                                        ...
            zeros(0,1),                                            ...
            zeros(0,1),                                            ...
            zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1),       ...
            zeros(0,1),                                            ... % safety cost
            zeroc,                                                 ...
            zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1),       ...
            zeros(0,1),                                            ... % constraint norm
            zeros(0,1), zeros(0,1),                                ...
            zeros(0,1), zeros(0,1),                                ...
            zeros(0,1), zeros(0,1),                                ...
            cell(0,1),                                             ...
            cell(0,1),                                             ... % per-obstacle distances
            'VariableNames', {                                     ...
                'iter',                                            ...
                'fval',                                            ...
                'length', 'curvature', 'time', 'jerk',            ...
                'safety',                                          ...
                'all_obs',                                         ...
                'g_obs',  'g_vel',     'g_acc', 'g_curv',         ...
                'g_norm',                                          ...
                'mean_speed', 'max_speed',                         ...
                'mean_acc',   'max_acc',                           ...
                'mean_curv',  'max_curv',                          ...
                'paths',                                           ...
                'obs_dists'                                        ...
            });
        disp(Constraint(x, Par))
    end

    if strcmp(state, 'iter')

        % ── Compute path and kinematics ───────────────────────────────────
        X = x(:) .* (ub - lb) + lb;
        P = bernstein_path(X, Par);
        [v, a, ~, kappa, ~, ~] = kinematics(P, Par);

        speed    = vecnorm(v, 2, 2);
        acc_norm = vecnorm(a, 2, 2);

        % ── Constraints ───────────────────────────────────────────────────
        [g, ~] = Constraint(x, Par);

        if iscell(Par.obs)
            n_obs = length(Par.obs);
        else
            n_obs = size(Par.obs, 1);
        end
        n_v = length(speed);
        n_a = length(acc_norm);

        g_obs  = max(g(1:n_obs));
        gg_obs = g(1:n_obs);
        g_vel  = max(g(n_obs+1         : n_obs+n_v));
        g_acc  = max(g(n_obs+n_v+1     : n_obs+n_v+n_a));
        g_curv = max(g(n_obs+n_v+n_a+1 : end));

        % Constraint vector norm (violated part only)
        g_norm = norm(max(g, 0));

        % ── Per-obstacle distances ────────────────────────────────────────
        [~, d_all] = obstacle_distance(P, Par);
        if numel(d_all) == n_obs
            obs_dists = d_all(:);
        else
            % Reconstruct from constraints: g = d_safe/d - 1  =>  d = d_safe/(g+1)
            obs_dists = Par.LengthReference * Par.d_safe ./ (gg_obs(:) + 1 + 1e-12);
        end

        % ── Cost components ───────────────────────────────────────────────
        if ~isempty(details)
            c_len    = details.length;
            c_curv   = details.curvature;
            c_time   = details.time;
            c_jerk   = details.jerk;
            c_safety = details.safety;
        else
            c_len = NaN; c_curv = NaN;
            c_time = NaN; c_jerk = NaN;
            c_safety = NaN;
        end

        % ── Append new row to table ───────────────────────────────────────
        new_row = table(                         ...
            optimValues.iteration,               ...
            optimValues.fval,                    ...
            c_len, c_curv, c_time, c_jerk,       ...
            c_safety,                            ...
            {gg_obs},                            ...
            g_obs, g_vel,  g_acc,  g_curv,       ...
            g_norm,                              ...
            mean(speed),   max(speed),            ...
            mean(acc_norm), max(acc_norm),        ...
            mean(abs(kappa)), max(abs(kappa)),    ...
            {P},                                 ...
            {obs_dists},                         ...
            'VariableNames', {                   ...
                'iter',                          ...
                'fval',                          ...
                'length', 'curvature', 'time', 'jerk', ...
                'safety',                        ...
                'all_obs',                       ...
                'g_obs',  'g_vel',     'g_acc', 'g_curv', ...
                'g_norm',                        ...
                'mean_speed', 'max_speed',       ...
                'mean_acc',   'max_acc',         ...
                'mean_curv',  'max_curv',        ...
                'paths',                         ...
                'obs_dists'                      ...
            });

        var_history = [var_history; new_row];

        % ── Build figure ──────────────────────────────────────────────────
        set(0, 'CurrentFigure', hfig);
        clf(hfig);

        % Layout: 4 rows x 6 cols  (tile numbering is row-major in tiledlayout)
        %   row1: 1  2  3  | 4  5  6
        %   row2: 7  8  9  | 10 11 12
        %   row3: 13 14 15 | 16 17 18
        %   row4: 19 20 21 | 22 23 24
        %
        % Map           : nexttile(1,[4,3])  → cols 1-3, all rows
        % Cost plots    : tiles 4,5,10,11   → cols 4-5, rows 1-2  (length,curvature,time,jerk)
        % Constr. plots : tiles 16,17,22,23 → cols 4-5, rows 3-4  (g_obs,g_vel,g_acc,g_curv)
        % NEW col 6     : tiles 6,12,18,24  → col 6,    rows 1-4  (fval,safety,g_norm,obs_dist)

        tl = tiledlayout(4, 6, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf('Degree: %d | Iter: %d | Cost: %.4f', ...
              n+1, optimValues.iteration, optimValues.fval), 'FontSize', 12);

        iters = var_history.iter;

        % ══════════════════════════════════════════════════════════════════
        % LEFT [4×3]: Trajectory map
        % ══════════════════════════════════════════════════════════════════
        ax_traj = nexttile(tl, 1, [4, 3]);
        hold(ax_traj, 'on'); grid(ax_traj, 'on');

        % Obstacles + buffer zones
        if iscell(Par.obs)
            buf = Par.d_safe;
            for i = 1:length(Par.obs)
                c = Par.obs{i};
                if size(c, 1) < 3, continue; end
                fill(ax_traj, c(:,1), c(:,2), [0.85 0.85 0.85], ...
                     'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.5);
                % Buffer zone: radially expand each vertex from centroid
                cx    = mean(c(:,1));  cy = mean(c(:,2));
                dc    = c - [cx cy];
                r     = vecnorm(dc, 2, 2);
                c_buf = c + dc ./ (r + 1e-12) * buf;
                plot(ax_traj, [c_buf(:,1); c_buf(1,1)], [c_buf(:,2); c_buf(1,2)], ...
                     '--', 'Color', [0.90 0.50 0.10], 'LineWidth', 1.0);
            end
        else
            theta = linspace(0, 2*pi, 60);
            buf   = Par.d_safe * Par.LengthReference;
            for i = 1:size(Par.obs, 1)
                ox = Par.obs(i,1) + Par.obs(i,3) * cos(theta);
                oy = Par.obs(i,2) + Par.obs(i,3) * sin(theta);
                fill(ax_traj, ox, oy, [0.85 0.85 0.85], ...
                     'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.5);
                bx = Par.obs(i,1) + (Par.obs(i,3) + buf) * cos(theta);
                by = Par.obs(i,2) + (Par.obs(i,3) + buf) * sin(theta);
                plot(ax_traj, bx, by, '--', ...
                     'Color', [0.90 0.50 0.10], 'LineWidth', 1.0);
            end
        end

        % All past paths as very thin grey ghost trails
        n_hist = height(var_history);
        for i = 1 : n_hist - 1
            alpha = 0.04 + 0.12 * (i / n_hist);
            Pg = var_history.paths{i};
            plot(ax_traj, Pg(:,1), Pg(:,2), '-', ...
                 'Color', [0.55 0.55 0.55 alpha], 'LineWidth', 0.5);
        end

        % Initial path on top of ghosts (dashed blue)
        if n_hist >= 1
            P_init = var_history.paths{1};
            plot(ax_traj, P_init(:,1), P_init(:,2), 'b--', 'LineWidth', 1.2);
        end

        % Current path
        plot(ax_traj, P(:,1), P(:,2), 'r-', 'LineWidth', 2);

        % Start / end markers
        plot(ax_traj, [Par.A(1), Par.B(1)], [Par.A(2), Par.B(2)], ...
             'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);

        axis(ax_traj, 'equal');
        xlabel(ax_traj, 'x [m]'); ylabel(ax_traj, 'y [m]');
        title(ax_traj, 'Trajectory');

        % ══════════════════════════════════════════════════════════════════
        % MIDDLE (cols 4-5): original 8 cost + constraint plots
        % ══════════════════════════════════════════════════════════════════

        % Cost histories — rows 1-2, cols 4-5  (tiles 4, 5, 10, 11)
        % All cost plots use blue tones for visual consistency
        cost_cfg = {
            4,  'length',    'Length',    [0.08 0.40 0.75], 'o';
            5,  'curvature', 'Curvature', [0.20 0.60 0.90], 's';
            10, 'time',      'Time',      [0.45 0.75 0.95], '^';
            11, 'jerk',      'Jerk',      [0.10 0.25 0.55], 'd';
        };

        for k = 1:size(cost_cfg, 1)
            ax = nexttile(tl, cost_cfg{k,1});
            hold(ax, 'on'); grid(ax, 'on');
            vals = var_history.(cost_cfg{k,2});
            if ~isempty(iters)
                plot(ax, iters, vals, ['-' cost_cfg{k,5}], ...
                     'Color', cost_cfg{k,4}, 'LineWidth', 1.5, 'MarkerSize', 3);
                xlim(ax, [0, max(iters)+1]);
            end
            title(ax,  cost_cfg{k,3}, 'FontSize', 9);
            xlabel(ax, 'Iter',        'FontSize', 8);
            ylabel(ax, 'Value',       'FontSize', 8);
            ax.FontSize = 8;
        end

        % Constraint histories — rows 3-4, cols 4-5  (tiles 16, 17, 22, 23)
        % All constraint plots use red tones for visual consistency
        con_cfg = {
            17, 'g_obs',  'Obstacle',  [0.80 0.10 0.10];
            22, 'g_vel',  'Velocity',  [0.90 0.35 0.20];
            23, 'g_acc',  'Accel.',    [0.70 0.15 0.30];
            24, 'g_curv', 'Curvature', [0.95 0.55 0.40];
        };

        for k = 1:size(con_cfg, 1)
            ax = nexttile(tl, con_cfg{k,1});
            hold(ax, 'on'); grid(ax, 'on');
            vals = var_history.(con_cfg{k,2});
            if ~isempty(iters)
                plot(ax, iters, vals, '-o', ...
                     'Color', con_cfg{k,4}, 'LineWidth', 1.5, 'MarkerSize', 3);
                last_val  = vals(end);
                dot_color = [0.20 0.75 0.30];
                if last_val > 0, dot_color = [0.90 0.20 0.20]; end
                plot(ax, iters(end), last_val, 'o', ...
                     'Color', dot_color, 'MarkerFaceColor', dot_color, 'MarkerSize', 5);
                xlim(ax, [0, max(iters)+1]);
            end
            yline(ax, 0, 'k--', 'LineWidth', 1.0);
            title(ax,  con_cfg{k,3}, 'FontSize', 9);
            xlabel(ax, 'Iter',       'FontSize', 8);
            ylabel(ax, 'max(g)',     'FontSize', 8);
            ax.FontSize = 8;
        end

        % ── Tile 6: Total cost f(x) ───────────────────────────────────────
        ax = nexttile(tl, 6);
        hold(ax, 'on'); grid(ax, 'on');
        plot(ax, iters, var_history.fval, '-o', ...
             'Color', [0.05 0.20 0.50], 'LineWidth', 1.5, 'MarkerSize', 3);
        if ~isempty(iters), xlim(ax, [0, max(iters)+1]); end
        title(ax,  'Total Cost', 'FontSize', 9);
        xlabel(ax, 'Iter',       'FontSize', 8);
        ylabel(ax, 'f(x)',       'FontSize', 8);
        ax.FontSize = 8;

        % ── Tile 12: Safety cost ──────────────────────────────────────────
        ax = nexttile(tl, 12);
        hold(ax, 'on'); grid(ax, 'on');
        plot(ax, iters, var_history.safety, '-o', ...
             'Color', [0.30 0.55 0.85], 'LineWidth', 1.5, 'MarkerSize', 3);
        if ~isempty(iters), xlim(ax, [0, max(iters)+1]); end
        title(ax,  'Safety Cost',  'FontSize', 9);
        xlabel(ax, 'Iter',         'FontSize', 8);
        ylabel(ax, 'c_{safety}',   'FontSize', 8);
        ax.FontSize = 8;

        % ── Tile 18: Constraint norm ||max(g,0)|| ────────────────────────
        ax = nexttile(tl, 16);
        hold(ax, 'on'); grid(ax, 'on');
        % g_norm is a scalar per iteration stored as a regular table column
        g_norm_hist = var_history.g_norm;
        plot(ax, iters, g_norm_hist, '-o', ...
             'Color', [0.60 0.05 0.15], 'LineWidth', 1.5, 'MarkerSize', 3);
        if ~isempty(iters)
            last_val  = g_norm_hist(end);
            dot_color = [0.20 0.75 0.30];
            if last_val > 0, dot_color = [0.90 0.20 0.20]; end
            plot(ax, iters(end), last_val, 'o', ...
                 'Color', dot_color, 'MarkerFaceColor', dot_color, 'MarkerSize', 5);
            xlim(ax, [0, max(iters)+1]);
        end
        yline(ax, 0, 'k--', 'LineWidth', 1.0);
        title(ax,  'Constraints norm',  'FontSize', 9);
        xlabel(ax, 'Iter',          'FontSize', 8);
        ylabel(ax, 'Constr. norm',  'FontSize', 8);
        ax.FontSize = 8;

        % ── Tile 24: Per-obstacle distance history (one line per obstacle) ──
        % Each obstacle gets its own line showing how its minimum clearance
        % evolves over iterations — analogous to the scalar constraint plots
        % but with n_obs traces stacked on the same axes.
        ax = nexttile(tl, 18);
        hold(ax, 'on'); grid(ax, 'on');
        n_hist = height(var_history);
        if n_hist > 0
            % Build matrix: rows = iterations, cols = obstacles
            n_obs_plot = numel(var_history.obs_dists{1});
            dist_mat   = NaN(n_hist, n_obs_plot);
            for ii = 1:n_hist
                d_row = var_history.obs_dists{ii};
                if numel(d_row) == n_obs_plot
                    dist_mat(ii, :) = d_row(:)';
                end
            end
            % Color map: blue shades, one per obstacle
            cmap = parula(n_obs_plot);
            for oi = 1:n_obs_plot
                plot(ax, iters, dist_mat(:, oi), '-', ...
                     'Color', cmap(oi,:), 'LineWidth', 1.0);
            end
            % Safe-distance reference line
            yline(ax, Par.d_safe, 'r--', 'LineWidth', 1.2);
            xlim(ax, [0, max(iters)+1]);
        end
        xlabel(ax, 'Iter',           'FontSize', 8);
        ylabel(ax, 'Min dist [m]',   'FontSize', 8);
        title(ax,  'Obs. clearances', 'FontSize', 9);
        ax.FontSize = 8;

        drawnow limitrate;
    end

    if strcmp(state, 'done')
        fprintf('Optimization complete. Final cost: %.4f\n', optimValues.fval);

        % Console summary (excluding cell columns)
        skip = ismember(var_history.Properties.VariableNames, {'paths','obs_dists','all_obs'});
        disp(var_history(:, ~skip));
    end
end