% MATLAB Script to visualize the 3D Cost Surface of a Bernstein Path
% This script iterates through values of alpha1 and alpha2 and plots the result

function generate_cost_surface()
    % --- 1. Parameters for the Grid Search ---
    resolution = 50;                  % Number of points per axis (50x50 = 2500 samples)
    limit = 20;                       % Search range for alphas [-20, 20]
    alpha_vals = linspace(-limit, limit, resolution);
    
    % Create the coordinate matrices for the parameter space
    [A1, A2] = meshgrid(alpha_vals, alpha_vals);
    
    % Initialize the matrix to store the cost results
    Z_cost = zeros(size(A1));
    
    fprintf('Computing cost surface for %d configurations...\n', resolution^2);

    % --- 2. Surface Computation ---
    % We loop through every combination of alpha1 and alpha2
    for i = 1:resolution
        for j = 1:resolution
            % Execute the path function provided by the user
            % Expected: function cost = bernstein_path(alpha1, alpha2)
            Z_cost(i,j) = bernstein_path(A1(i,j), A2(i,j));
        end
    end

    % --- 3. 3D Visualization ---
    fig = figure('Color', 'w', 'Name', 'Bernstein Path Cost Landscape');
    
    % Render the 3D surface
    surf_plot = surf(A1, A2, Z_cost);
    
    % Aesthetics and Rendering
    shading interp;          % Interpolate colors for a smooth look
    colormap(jet);           % Classic heat map (Blue = Low Cost, Red = High Cost)
    colorbar;                % Show the cost scale
    camlight;                % Add a light source to highlight the surface curvature
    
    % Labels and Annotations
    xlabel('\alpha_1 (Initial bending)');
    ylabel('\alpha_2 (Final bending)');
    zlabel('Total Cost J');
    title('3D Trajectory Energy Landscape');
    
    % --- 4. Finding the Global Minimum ---
    % Locate the absolute lowest cost on the surface
    [min_val, min_idx] = min(Z_cost(:));
    [row, col] = ind2sub(size(Z_cost), min_idx);
    
    hold on;
    % Mark the optimal parameter pair with a red sphere
    plot3(A1(row,col), A2(row,col), min_val, 'ro', ...
          'MarkerSize', 12, 'MarkerFaceColor', 'r');
    
    text(A1(row,col), A2(row,col), min_val + (0.05 * max(Z_cost(:))), ...
         'Global Optimum', 'Color', 'r', 'FontWeight', 'bold');
    
    % Set the camera angle for best perspective
    view(-45, 30);
    grid on;
end