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
            zeroc,                                              ...
            zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1),       ...
            zeros(0,1), zeros(0,1),                                ...
            zeros(0,1), zeros(0,1),                                ...
            zeros(0,1), zeros(0,1),                                ...
            cell(0,1),                                             ...
            'VariableNames', {                                     ...
                'iter',                                            ...
                'fval',                                            ...
                'length', 'curvature', 'time', 'jerk',            ...
                'all_obs',                                        ...
                'g_obs',  'g_vel',     'g_acc', 'g_curv',         ...
                'mean_speed', 'max_speed',                         ...
                'mean_acc',   'max_acc',                           ...
                'mean_curv',  'max_curv',                          ...
                'paths'                                            ...
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

        % ── Cost components ───────────────────────────────────────────────
        if ~isempty(details)
            c_len  = details.length;
            c_curv = details.curvature;
            c_time = details.time;
            c_jerk = details.jerk;
        else
            c_len  = NaN; c_curv = NaN;
            c_time = NaN; c_jerk = NaN;
        end

        % ── Append new row to table ───────────────────────────────────────
        new_row = table(                         ...
            optimValues.iteration,               ...
            optimValues.fval,                    ...
            c_len, c_curv, c_time, c_jerk,       ...
            {gg_obs},                              ... 
            g_obs, g_vel,  g_acc,  g_curv,       ...
            mean(speed),   max(speed),            ...
            mean(acc_norm), max(acc_norm),        ...
            mean(abs(kappa)), max(abs(kappa)),    ...
            {P},                                 ...   % cell-wrap for table
            'VariableNames', {                   ...
                'iter',                          ...
                'fval',                          ...
                'length', 'curvature', 'time', 'jerk',     ...
                'all_obs'                                  ...
                'g_obs',  'g_vel',     'g_acc', 'g_curv',  ...
                'mean_speed', 'max_speed',       ...
                'mean_acc',   'max_acc',         ...
                'mean_curv',  'max_curv',        ...
                'paths'                          ...
            });

        var_history = [var_history; new_row];

        % ── Build figure ──────────────────────────────────────────────────
        set(0, 'CurrentFigure', hfig);
        clf(hfig);

        tl = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf('Degree: %d | Iter: %d | Cost: %.4f', ...
              n+1, optimValues.iteration, optimValues.fval), 'FontSize', 12);

        iters = var_history.iter;

        % ── LEFT HALF [4x2]: Trajectory ──────────────────────────────────
        ax_traj = nexttile(tl, 1, [4, 2]);
        hold(ax_traj, 'on'); grid(ax_traj, 'on');

        if iscell(Par.obs)
            for i = 1:length(Par.obs)
                c = Par.obs{i};
                if size(c, 1) < 3, continue; end
                fill(ax_traj, c(:,1), c(:,2), [0.85 0.85 0.85], ...
                     'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.5);
            end
        else
            theta = linspace(0, 2*pi, 60);
            for i = 1:size(Par.obs, 1)
                ox = Par.obs(i,1) + Par.obs(i,3) * cos(theta);
                oy = Par.obs(i,2) + Par.obs(i,3) * sin(theta);
                fill(ax_traj, ox, oy, [0.85 0.85 0.85], ...
                     'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.5);
            end
        end

        % Ghost trails: last 5 paths faded
        n_hist = height(var_history);
        for i = max(1, n_hist-4) : n_hist-1
            alpha = 0.10 + 0.15 * (i - max(1, n_hist-4));
            Pg = var_history.paths{i};
            plot(ax_traj, Pg(:,1), Pg(:,2), '-', ...
                 'Color', [0.6 0.6 0.6 alpha], 'LineWidth', 1);
        end

        plot(ax_traj, P(:,1), P(:,2), 'r-', 'LineWidth', 2);
        plot(ax_traj, [Par.A(1), Par.B(1)], [Par.A(2), Par.B(2)], ...
             'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
        axis(ax_traj, 'equal');
        xlabel(ax_traj, 'x [m]'); ylabel(ax_traj, 'y [m]');
        title(ax_traj, 'Trajectory');

        % ── UP-RIGHT [2x2]: Cost histories ───────────────────────────────
        cost_cfg = {
            3,  'length',    'Length',    [0.85 0.33 0.10], 'o';
            4,  'curvature', 'Curvature', [0.47 0.67 0.19], 's';
            7,  'time',      'Time',      [0.30 0.60 0.90], '^';
            8,  'jerk',      'Jerk',      [0.60 0.20 0.80], 'd';
        };

        for k = 1:size(cost_cfg, 1)
            ax = nexttile(tl, cost_cfg{k,1});
            hold(ax, 'on'); grid(ax, 'on');
            vals = var_history.(cost_cfg{k,2});   % table column access
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

        % ── DOWN-RIGHT [2x2]: Constraint histories ────────────────────────
        con_cfg = {
            11, 'g_obs',  'Obstacle';
            12, 'g_vel',  'Velocity';
            15, 'g_acc',  'Accel.';
            16, 'g_curv', 'Curvature';
        };

        for k = 1:size(con_cfg, 1)
            ax = nexttile(tl, con_cfg{k,1});
            hold(ax, 'on'); grid(ax, 'on');
            vals = var_history.(con_cfg{k,2});    % table column access
            if ~isempty(iters)
                plot(ax, iters, vals, '-o', ...
                     'Color', [0.30 0.60 0.90], 'LineWidth', 1.5, 'MarkerSize', 3);
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

        drawnow limitrate;
    end

    if strcmp(state, 'done')
        fprintf('Optimization complete. Final cost: %.4f\n', optimValues.fval);

        % Console summary (excluding paths column)
        disp(var_history(:, ~strcmp(var_history.Properties.VariableNames, 'paths')));
    end
end