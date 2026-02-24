% MATLAB Interactive Bernstein Path Planner
% Parameters: A, B (Start/End), alpha1, alpha2 (Shape variables)

function interactive_bernstein()
    % --- Initial Configuration ---
    A = [0; 0];      % Start Point
    B = [10; 0];     % End Point
    
    % Calculate the vector AB and its perpendicular normal vector
    dir_vec = B - A;
    L = norm(dir_vec);
    % Normal vector (rotate 90 degrees)
    n_vec = [-dir_vec(2); dir_vec(1)] / L;
    
    % Time vector for the curve (0 to 1)
    t = linspace(0, 1, 100);
    
    % Precompute Bernstein Basis Functions (Cubic)
    % B1,3 = 3*t*(1-t)^2
    % B2,3 = 3*t^2*(1-t)
    b1 = 3 .* t .* (1-t).^2;
    b2 = 3 .* t.^2 .* (1-t);
    
    % --- Setup UI Figure ---
    fig = figure('Name', 'Interactive Bernstein Path', 'Position', [100, 100, 800, 600]);
    ax = axes('Parent', fig, 'Position', [0.1, 0.3, 0.8, 0.6]);
    grid(ax, 'on'); hold(ax, 'on'); axis(ax, 'equal');
    xlabel('X Position'); ylabel('Y Position');
    title('Adjust \alpha_1 and \alpha_2 to deform the path');
    xlim([-2, 12]); ylim([-5, 5]);

    % Initial Plot handles
    path_plot = plot(ax, 0, 0, 'b', 'LineWidth', 2);
    plot(ax, A(1), A(2), 'ro', 'MarkerFaceColor', 'r'); % Start
    plot(ax, B(1), B(2), 'go', 'MarkerFaceColor', 'g'); % End
    
    % --- UI Controls (Sliders) ---
    uicontrol('Style', 'text', 'Position', [150, 80, 100, 20], 'String', 'Alpha 1');
    s1 = uicontrol('Style', 'slider', 'Min', -10, 'Max', 10, 'Value', 0, ...
                   'Position', [250, 80, 300, 20], 'Callback', @update_plot);
               
    uicontrol('Style', 'text', 'Position', [150, 40, 100, 20], 'String', 'Alpha 2');
    s2 = uicontrol('Style', 'slider', 'Min', -10, 'Max', 10, 'Value', 0, ...
                   'Position', [250, 40, 300, 20], 'Callback', @update_plot);

    % Call update once to show initial line
    update_plot();

    % --- Nested Update Function ---
    function update_plot(~, ~)
        % Get values from sliders
        val1 = get(s1, 'Value');
        val2 = get(s2, 'Value');
        
        % Calculate Baseline (Linear interpolation between A and B)
        baseline = A + t .* dir_vec;
        
        % Calculate Deviation (Perpendicular component)
        % deviation = n_vec * (alpha1*b1 + alpha2*b2)
        deviation = n_vec * (val1 .* b1 + val2 .* b2);
        
        % Resulting Curve
        P = baseline + deviation;
        
        % Update Plot
        set(path_plot, 'XData', P(1,:), 'YData', P(2,:));
    end
end