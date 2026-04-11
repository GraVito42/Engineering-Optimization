function Plotinitialguess(Par, x0s, idx)
% Plotinitialguess  Visualise all sampled initial paths and the chosen one.
%
%   PlotInitialGuesses(Par, x0s, idx, ub, lb)
%
%   Inputs
%   ------
%   Par  : parameter struct (same as used by cost_function / Constraint)
%   x0s  : (N x n_vars) matrix – every sampled initial guess in [0,1] space
%   idx  : scalar row index into x0s that was selected as the best guess
%   ub   : upper bounds vector  (n_vars x 1)
%   lb   : lower bounds vector  (n_vars x 1)
%
%   The function draws:
%     • All N candidate paths in semi-transparent blue
%     • The selected path (x0s(idx,:)) in solid red
%     • Start / end waypoints as filled circles
%     • Obstacles (both the cell-polygon and the matrix-circle formats)

    global ub lb

    n_candidates = size(x0s, 1);

    figure('Name', 'Initial Guess Survey', 'Color', 'w', ...
           'Units', 'normalized', 'Position', [0.05 0.1 0.55 0.75]);

    tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('Initial Guess Survey  |  %d candidates  |  selected idx = %d', ...
          n_candidates, idx), 'FontSize', 13, 'FontWeight', 'bold');

    ax = nexttile(tl);
    hold(ax, 'on');
    grid(ax, 'on');
    axis(ax, 'equal');
    xlabel(ax, 'x [m]', 'FontSize', 11);
    ylabel(ax, 'y [m]', 'FontSize', 11);
    title(ax, 'Trajectory Overview', 'FontSize', 11);

    % ------------------------------------------------------------------
    % 1.  Obstacles
    % ------------------------------------------------------------------
    if iscell(Par.obs)
        % Polygon contours from island_detection()
        for i = 1:length(Par.obs)
            c = Par.obs{i};
            if size(c, 1) < 3, continue; end
            fill(ax, c(:,1), c(:,2), [0.82 0.82 0.82], ...
                 'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.55, ...
                 'HandleVisibility', 'off');
        end
    else
        % Circular obstacles  [cx, cy, r]
        theta = linspace(0, 2*pi, 60);
        for i = 1:size(Par.obs, 1)
            ox = Par.obs(i,1) + Par.obs(i,3) * cos(theta);
            oy = Par.obs(i,2) + Par.obs(i,3) * sin(theta);
            fill(ax, ox, oy, [0.82 0.82 0.82], ...
                 'EdgeColor', [0.4 0.4 0.4], 'FaceAlpha', 0.55, ...
                 'HandleVisibility', 'off');
        end
    end

    % ------------------------------------------------------------------
    % 2.  All candidate paths  (blue, semi-transparent)
    % ------------------------------------------------------------------
    h_all = [];
    for k = 1:n_candidates
        X = x0s(k,:) .* (ub - lb)' + lb';   % de-normalise to real coords
        P = bernstein_path(X, Par);
        h = plot(ax, P(:,1), P(:,2), '-', ...
                 'Color', [0.20 0.45 0.85 0.18], ...   % RGBA – low alpha
                 'LineWidth', 0.8, ...
                 'HandleVisibility', 'off');
        if k == 1
            h_all = h;                                  % keep one for legend
            set(h_all, 'HandleVisibility', 'on');
        end
    end

    % ------------------------------------------------------------------
    % 3.  Selected path  (red, solid, on top)
    % ------------------------------------------------------------------
    X_best = x0s(idx,:) .* (ub - lb)' + lb';
    P_best = bernstein_path(X_best, Par);
    h_best = plot(ax, P_best(:,1), P_best(:,2), 'r-', ...
                  'LineWidth', 2.5, 'DisplayName', 'Selected path');

    % ------------------------------------------------------------------
    % 4.  Start / end markers
    % ------------------------------------------------------------------
    plot(ax, Par.A(1), Par.A(2), 'ko', ...
         'MarkerFaceColor', [0.10 0.80 0.10], 'MarkerSize', 10, ...
         'DisplayName', 'Start');
    plot(ax, Par.B(1), Par.B(2), 'ko', ...
         'MarkerFaceColor', [0.90 0.20 0.20], 'MarkerSize', 10, ...
         'DisplayName', 'Goal');

    % ------------------------------------------------------------------
    % 5.  Legend
    % ------------------------------------------------------------------
    set(h_all, 'Color', [0.20 0.45 0.85 0.6], 'LineWidth', 1.2, ...
               'DisplayName', sprintf('Candidates (%d)', n_candidates));
    legend(ax, 'Location', 'best', 'FontSize', 10);

    hold(ax, 'off');
    drawnow;
end