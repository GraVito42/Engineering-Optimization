function stop = TrajectoryPlotter(x, optimValues, state)
    stop = false;
    global Par 
    global details
    global ub lb
    n = length(x);
    
    % Figure 1: 2D Trajectory Animation
    figure(1); 
    if strcmp(state, 'iter')
        clf;
        tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf('Degree: %d | Iter: %d | Cost: %.4f', ...
              n+1, optimValues.iteration, optimValues.fval), 'FontSize', 12);

        % ── Trajectory panel ──────────────────────────────────────────────
        ax_traj = nexttile(tl, 1, [2, 2]);
        hold(ax_traj, 'on'); grid(ax_traj, 'on');

        % ── Obstacle rendering: supports both formats ─────────────────────
        if iscell(Par.obs)
            % NEW: generic polygon contours from island_detection()
            for i = 1:length(Par.obs)
                c = Par.obs{i};
                if size(c, 1) < 3, continue; end
                fill(ax_traj, c(:,1), c(:,2), [0.85 0.85 0.85], ...
                     'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.5);
            end
        else
            % OLD: circular obstacles [cx, cy, r]
            theta = linspace(0, 2*pi, 60);
            for i = 1:size(Par.obs, 1)
                ox = Par.obs(i,1) + Par.obs(i,3) * cos(theta);
                oy = Par.obs(i,2) + Par.obs(i,3) * sin(theta);
                fill(ax_traj, ox, oy, [0.85 0.85 0.85], ...
                     'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.5);
            end
        end
        X = x(:) .* (ub - lb) + lb; % Scale from [0,1] to actual bounds
        P = bernstein_path(X, Par);
        plot(ax_traj, P(:,1), P(:,2), 'r-', 'LineWidth', 2);
        plot(ax_traj, [Par.A(1), Par.B(1)], [Par.A(2), Par.B(2)], ...
             'ko', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
        axis(ax_traj, 'equal');
        xlabel(ax_traj, 'x [m]'); ylabel(ax_traj, 'y [m]');
        title(ax_traj, 'Trajectory');

        % ── Cost bar chart ────────────────────────────────────────────────
        ax_cost = nexttile(tl, 3);
        if ~isempty(details)
            cost_names = {'Length', 'Curv.', 'Safety', 'Time', 'Jerk'};
            cost_vals  = [details.length, details.curvature, details.safety, ...
                          details.time,   details.jerk];
            b = bar(ax_cost, cost_vals', 'grouped');
            b(1).FaceColor = [0.30 0.60 0.90];
            set(ax_cost, 'XTickLabel', cost_names, 'FontSize', 9);
            legend(ax_cost, 'Raw', 'Location', 'northwest', 'FontSize', 8);
            title(ax_cost, 'Cost Components');
            ylabel(ax_cost, 'Value');
            grid(ax_cost, 'on');
        end

        % ── Constraint bar chart ──────────────────────────────────────────
        ax_con = nexttile(tl, 6);
        [g, ~] = Constraint(x, Par);

        % Obstacle constraint count depends on Par.obs format
        if iscell(Par.obs)
            n_obs = length(Par.obs);
        else
            n_obs = size(Par.obs, 1);
        end

        [v, a, ~, k, ~, ~] = kinematics(P, Par);
        n_v = length(vecnorm(v, 2, 2));
        n_a = length(vecnorm(a, 2, 2));
        n_k = length(k);

        g_obs = max(g(1:n_obs));
        g_v   = max(g(n_obs+1         : n_obs+n_v));
        g_a   = max(g(n_obs+n_v+1     : n_obs+n_v+n_a));
        g_k   = max(g(n_obs+n_v+n_a+1 : end));

        g_summary = [g_obs, g_v, g_a, g_k];
        labels_s  = {'Obstacle', 'Velocity', 'Accel.', 'Curvature'};
        colors    = repmat([0.20 0.75 0.30], 4, 1);
        colors(g_summary > 0, :) = repmat([0.90 0.20 0.20], sum(g_summary > 0), 1);

        bh = bar(ax_con, g_summary, 'FaceColor', 'flat');
        bh.CData = colors;
        hold(ax_con, 'on');
        yline(ax_con, 0, 'k--', 'LineWidth', 1.2);
        set(ax_con, 'XTick', 1:4, 'XTickLabel', labels_s, 'FontSize', 9);
        title(ax_con, 'Constraints');
        ylabel(ax_con, 'max(g)');
        grid(ax_con, 'on');

        drawnow limitrate;

        % ── 3D surface plot (n=2 only) ────────────────────────────────────
        if n == 2
            figure(2);
            if strcmp(state, 'init')
                Par.Surf = GenerateSurfaceData(Par);
                surf(Par.Surf.A1, Par.Surf.A2, Par.Surf.Z, ...
                     'EdgeColor', 'none', 'FaceAlpha', 0.6);
                hold on; colormap jet; view(-45, 30);
            elseif strcmp(state, 'iter')
                plot3(x(1), x(2), optimValues.fval, 'k.', 'MarkerSize', 5);
            end
        end
    end
    
    if strcmp(state, 'done')
        fprintf('Optimization complete. Final Cost: %.4f\n', optimValues.fval);
    end
end
