function stop = TrajectoryPlotter(x, optimValues, state)
    stop = false;
    global Par
    n = length(x);
    
    % Figure 1: 2D Trajectory Animation
    figure(1); 
    if strcmp(state, 'iter')
        hold on; grid on; cla;
        
        % Plot Obstacles
        theta = linspace(0, 2*pi, 30);
        for i = 1:size(Par.obs, 1)
            ox = Par.obs(i,1) + Par.obs(i,3) * cos(theta);
            oy = Par.obs(i,2) + Par.obs(i,3) * sin(theta);
            fill(ox, oy, [0.8 0.8 0.8], 'EdgeColor', 'k', 'FaceAlpha', 0.4);
        end
        
        P = bernstein_path(x, Par);
        plot(P(:,1), P(:,2), 'r-', 'LineWidth', 1.5);
        plot([Par.A(1), Par.B(1)], [Par.A(2), Par.B(2)], 'ko', 'MarkerFaceColor', 'g');
        title(sprintf('Degree: %d | Iter: %d | Cost: %.4f', n+1, optimValues.iteration, optimValues.fval));
        axis equal; drawnow limitrate;
    end

    % Figure 2: 3D Surface Plot (Only for n=2)
    if n == 2
        figure(2);
        if strcmp(state, 'init')
            % Pre-compute Surface once (as discussed previously)
            Par.Surf = GenerateSurfaceData(Par); 
            surf(Par.Surf.A1, Par.Surf.A2, Par.Surf.Z, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
            hold on; colormap jet; view(-45, 30);
        elseif strcmp(state, 'iter')
            plot3(x(1), x(2), optimValues.fval, 'k.', 'MarkerSize', 5);
        end
    end
    
    if strcmp(state, 'done')
        fprintf('Optimization complete. Final Cost: %.4f\n', optimValues.fval);
    end
end