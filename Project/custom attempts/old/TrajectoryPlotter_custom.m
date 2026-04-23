function stop = TrajectoryPlotter_custom(x, optimValues, state)
% TrajectoryPlotter  OutputFcn for fmincon.
%
%   On 'init'  : draws all initial-guess candidate paths (blue, faded) as a
%                permanent background layer, then overlays the starting path.
%   On 'iter'  : updates only the live trajectory + cost/constraint bars
%                without wiping the candidate-path background.
%   On 'done'  : prints a summary to the command window.
%
%   Requires globals: Par, details, ub, lb, x0s, idx
%     x0s  – (N x n_vars) matrix of ALL sampled initial guesses in [0,1]
%     idx  – index of the chosen initial guess (used to highlight it)

    stop = false;
    global Par details ub lb x0s idx

    n = length(x);

    % ======================================================================
    %  STATE: init
    %  Build figure skeleton + draw candidate background once.
    % ======================================================================
    if strcmp(state, 'init')

        figure(1);
        clf;
        tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl, sprintf('Degree: %d | Iter: 0 | Cost: –', n+1), 'FontSize', 12);

        % ── Trajectory panel (spans left 2 columns) ───────────────────────
        ax_traj = nexttile(tl, 1, [2, 2]);
        hold(ax_traj, 'on');
        grid(ax_traj, 'on');
        axis(ax_traj, 'equal');
        xlabel(ax_traj, 'x [m]');
        ylabel(ax_traj, 'y [m]');
        title(ax_traj, 'Trajectory');

        % ── Obstacles ─────────────────────────────────────────────────────
        draw_obstacles(ax_traj, Par);

        % ── Candidate paths (permanent background) ────────────────────────
        if ~isempty(x0s)
            n_cand = size(x0s, 1);
            h_first = [];
            for k = 1:n_cand
                X = x0s(k,:) .* (ub - lb)' + lb';
                P = bernstein_path(X, Par);
                hc = plot(ax_traj, P(:,1), P(:,2), '-', ...
                          'Color',     [0.20 0.45 0.85 0.15], ...
                          'LineWidth', 0.7, ...
                          'Tag',       'candidate', ...
                          'HandleVisibility', 'off');
                if k == 1, h_first = hc; end
            end

            % Highlight the selected initial guess in orange
            X_sel = x0s(idx,:) .* (ub - lb)' + lb';
            P_sel = bernstein_path(X_sel, Par);
            plot(ax_traj, P_sel(:,1), P_sel(:,2), '-', ...
                 'Color',     [1.00 0.55 0.00 0.85], ...
                 'LineWidth', 1.8, ...
                 'Tag',       'selected_guess', ...
                 'DisplayName', sprintf('Selected x_0 (idx %d)', idx));

            % Give the blue group one legend entry
            if ~isempty(h_first)
                set(h_first, 'HandleVisibility', 'on', ...
                    'DisplayName', sprintf('Candidates (%d)', n_cand));
            end
        end

        % ── Start / goal markers ──────────────────────────────────────────
        plot(ax_traj, Par.A(1), Par.A(2), 'ko', ...
             'MarkerFaceColor', [0.10 0.80 0.10], 'MarkerSize', 9, ...
             'DisplayName', 'Start');
        plot(ax_traj, Par.B(1), Par.B(2), 'ko', ...
             'MarkerFaceColor', [0.90 0.20 0.20], 'MarkerSize', 9, ...
             'DisplayName', 'Goal');

        legend(ax_traj, 'Location', 'best', 'FontSize', 8);

        % ── Reserve cost & constraint axes (empty for now) ────────────────
        ax_cost = nexttile(tl, 3);
        title(ax_cost, 'Cost Components'); grid(ax_cost, 'on');

        ax_con = nexttile(tl, 6);
        title(ax_con, 'Constraints'); grid(ax_con, 'on');

        drawnow;
        return          % nothing more to do at init
    end

    % ======================================================================
    %  STATE: iter
    %  Refresh ONLY the live trajectory line + the two bar charts.
    %  The candidate background (Tag='candidate') is left untouched.
    % ======================================================================
    if strcmp(state, 'iter')

        figure(1);

        % Retrieve the tiled layout and its children
        fig   = gcf;
        tl    = fig.Children(end);     % tiledlayout is the last child added

        % Update figure-level title
        tl.Title.String = sprintf('Degree: %d | Iter: %d | Cost: %.4f', ...
                                   n+1, optimValues.iteration, optimValues.fval);

        % ── Trajectory axis: delete only the live-path line ───────────────
        ax_traj = nexttile(tl, 1);

        % Remove previous live trajectory (tagged 'live_path')
        delete(findobj(ax_traj, 'Tag', 'live_path'));

        % Draw updated trajectory
        X = x(:) .* (ub - lb) + lb;
        P = bernstein_path(X, Par);
        plot(ax_traj, P(:,1), P(:,2), 'r-', ...
             'LineWidth', 2.2, ...
             'Tag', 'live_path', ...
             'HandleVisibility', 'off');

        % ── Cost bar chart ────────────────────────────────────────────────
        ax_cost = nexttile(tl, 3);
        cla(ax_cost);
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
        cla(ax_con);

        [g, ~] = Constraint(x, Par);

        if iscell(Par.obs)
            n_obs = length(Par.obs);
        else
            n_obs = size(Par.obs, 1);
        end

        [v, a, ~, k_cur, ~, ~] = kinematics(P, Par);
        n_v = length(vecnorm(v, 2, 2));
        n_a = length(vecnorm(a, 2, 2));

        g_obs = max(g(1 : n_obs));
        g_v   = max(g(n_obs+1          : n_obs+n_v));
        g_a   = max(g(n_obs+n_v+1      : n_obs+n_v+n_a));
        g_k   = max(g(n_obs+n_v+n_a+1  : end));

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

        % ── Optional 3-D surface (n=2 only) ──────────────────────────────
        if n == 2
            figure(2);
            plot3(x(1), x(2), optimValues.fval, 'k.', 'MarkerSize', 5);
            drawnow limitrate;
        end
    end

    % ======================================================================
    %  STATE: done
    % ======================================================================
    if strcmp(state, 'done')
        fprintf('Optimization complete. Final cost: %.6f\n', optimValues.fval);
    end
end


% =========================================================================
%  Local helper – obstacle rendering (shared between init & iter)
% =========================================================================
function draw_obstacles(ax, Par)
    if iscell(Par.obs)
        for i = 1:length(Par.obs)
            c = Par.obs{i};
            if size(c, 1) < 3, continue; end
            fill(ax, c(:,1), c(:,2), [0.85 0.85 0.85], ...
                 'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.50, ...
                 'HandleVisibility', 'off');
        end
    else
        theta = linspace(0, 2*pi, 60);
        for i = 1:size(Par.obs, 1)
            ox = Par.obs(i,1) + Par.obs(i,3) * cos(theta);
            oy = Par.obs(i,2) + Par.obs(i,3) * sin(theta);
            fill(ax, ox, oy, [0.85 0.85 0.85], ...
                 'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.50, ...
                 'HandleVisibility', 'off');
        end
    end
end