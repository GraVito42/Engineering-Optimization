% MATLAB Dynamic Bernstein Path Planner
% Allows for N-number of shape variables

function dynamic_bernstein()
    % --- Configuration ---
    A = [0; 0]; 
    B = [10; 0];
    n_vars = 8; % <--- CHANGE THIS: Number of independent variables
    
    dir_vec = B - A;
    L = norm(dir_vec);
    n_vec = [-dir_vec(2); dir_vec(1)] / L; % Normal vector
    t = linspace(0, 1, 100)'; % Column vector
    
    % --- Generate Bernstein Basis Matrix ---
    % Degree of polynomial d = n_vars + 1
    d = n_vars + 1;
    M = zeros(length(t), n_vars);
    for i = 1:n_vars
        % Formula: nchoosek(d, i) * t^i * (1-t)^(d-i)
        coeff = nchoosek(d, i);
        M(:, i) = coeff .* (t.^i) .* ((1-t).^(d-i));
    end
    
    % --- Setup UI Figure ---
    fig = figure('Name', sprintf('Dynamic Bernstein with %d vars', n_vars));
    ax = axes('Parent', fig, 'Position', [0.1, 0.4, 0.8, 0.5]);
    grid on; hold on; axis equal;
    xlim([-2, 12]); ylim([-5, 5]);
    path_plot = plot(ax, 0, 0, 'b', 'LineWidth', 2);
    
    % --- Dynamic Slider Generation ---
    sliders = [];
    spacing = 0.25 / n_vars;
    for i = 1:n_vars
        y_pos = 50 + (i-1)*(300/n_vars);
        uicontrol('Style', 'text', 'Position', [50, y_pos, 80, 20], ...
                  'String', sprintf('Alpha %d', i));
        sliders(i) = uicontrol('Style', 'slider', 'Min', -15, 'Max', 15, 'Value', 0, ...
                               'Position', [140, y_pos, 500, 20], 'Callback', @update_plot);
    end

    update_plot();

    % --- Nested Update Function ---
    function update_plot(~, ~)
        % Collect all alpha values into a vector
        alphas = zeros(1, n_vars);
        for k = 1:n_vars
            alphas(k) = get(sliders(k), 'Value');
        end
        
        % Calculate total deviation: Matrix multiplication M * alphas'
        deviation_mag = M * alphas'; % Result is [100x1]
        
        % Global position = Straight line + Normal deviation
        baseline = A' + t * dir_vec';
        P = baseline + deviation_mag * n_vec';
        
        set(path_plot, 'XData', P(:,1), 'YData', P(:,2));
    end
end