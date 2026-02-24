function stop = TrajectoryPlotter(x, optimValues, state)
    stop = false;
    global Par
    persistent history surfData
    
    % --- 1. INITIALIZATION ---
    if strcmp(state, 'init')
        history = [];
        % Pre-compute Surface for 3D plot
        res = 500; 
        % It's safer to pull bounds directly from the workspace if possible,
        % otherwise, ensure these match your lb/ub in main.m
        a_vec = linspace(-15, 15, res); 
        [A1, A2] = meshgrid(a_vec, a_vec);
        Z = zeros(size(A1));
        
        % Pre-calculate the landscape (cost + feasibility "cuts")
        for i = 1:numel(A1)
            alpha_t = [A1(i), A2(i)];
            Z(i) = Objective(alpha_t, Par);
            % Use the actual constraint function to mask infeasible regions
            if any(Constraint(alpha_t, Par) > 0)
                Z(i) = NaN; 
            end
        end
        surfData.A1 = A1; 
        surfData.A2 = A2; 
        surfData.Z = Z;
    end

    % --- 2. ITERATION (Animation) ---
    if strcmp(state, 'iter')
        % Track the optimization path in design space
        history = [history; x(1), x(2), optimValues.fval];
        
        % Figure 100: Live 2D Trajectory
        figure(100); clf; hold on; grid on;
        theta = linspace(0, 2*pi, 20);
        for i = 1:size(Par.obs,1)
            fill(Par.obs(i,1)+Par.obs(i,3)*cos(theta), Par.obs(i,2)+Par.obs(i,3)*sin(theta), ...
                [0.8 0.8 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.4);
        end
        P = bernstein_path(x(1), x(2), Par);
        plot(P(:,1), P(:,2), 'r-', 'LineWidth', 2);
        plot([Par.A(1), Par.B(1)], [Par.A(2), Par.B(2)], 'go', 'MarkerFaceColor', 'g');
        title(sprintf('Trajectory Iter: %d | Cost: %.4f', optimValues.iteration, optimValues.fval));
        axis equal; drawnow limitrate;

        % Figure 200: Live 3D Design Space
        figure(200); clf; hold on;
        surf(surfData.A1, surfData.A2, surfData.Z, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        if ~isempty(history)
            plot3(history(:,1), history(:,2), history(:,3), 'r-', 'LineWidth', 1.5);
        end
        plot3(x(1), x(2), optimValues.fval, 'ko', 'MarkerFaceColor', 'y');
        colormap jet; view(-45, 30); grid on;
        title('Design Space (Cost Landscape)');
        xlabel('\alpha_1'); ylabel('\alpha_2'); zlabel('Cost');
    end

    % --- 3. FINAL STATE (Summary) ---
    if strcmp(state, 'done')
        % Create a high-quality summary window
        figure('Name', 'Optimization Summary', 'Color', 'w');
        
        % Subplot 1: Final Trajectory with Obstacles
        subplot(1,2,1); hold on; grid on;
        theta = linspace(0, 2*pi, 50);
        for i = 1:size(Par.obs, 1)
            % Draw obstacles in a darker shade for the final report
            fill(Par.obs(i,1)+Par.obs(i,3)*cos(theta), Par.obs(i,2)+Par.obs(i,3)*sin(theta), ...
                [0.3 0.3 0.3], 'EdgeColor', 'k', 'FaceAlpha', 0.7);
        end
        P_final = bernstein_path(x(1), x(2), Par);
        plot(P_final(:,1), P_final(:,2), 'b-', 'LineWidth', 3);
        plot([Par.A(1), Par.B(1)], [Par.A(2), Par.B(2)], 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'y');
        title('Final Optimal Path'); 
        xlabel('X [m]'); ylabel('Y [m]');
        axis equal;

        % Subplot 2: Design Space Convergence
        subplot(1,2,2); hold on; grid on;
        surf(surfData.A1, surfData.A2, surfData.Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
        if ~isempty(history)
            plot3(history(:,1), history(:,2), history(:,3), 'k-', 'LineWidth', 1.5);
        end
        % Mark the converged point with a green star
        plot3(x(1), x(2), optimValues.fval, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
        title('Convergence Path'); 
        xlabel('\alpha_1'); ylabel('\alpha_2');
        view(-45, 30); colormap jet;
    end
end