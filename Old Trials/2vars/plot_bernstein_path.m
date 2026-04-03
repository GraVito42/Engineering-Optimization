function plot_bernstein_path(alpha1, alpha2)
    % PLOT_BERNSTEIN_PATH Plots a path and its associated cost
    % Inputs: alpha1, alpha2 (Cubic Bernstein coefficients)
    
    % --- Configuration ---
    A = [0; 0];      % Start Point
    B = [10; 0];     % End Point
    dt = 0.01;       % Time step
    t = (0:dt:1)';   % Time vector (column)
    
    % Geometry setup
    dir_vec = B - A;
    L = norm(dir_vec);
    n_vec = [-dir_vec(2); dir_vec(1)] / L; % Normal vector
    
    % --- Basis Computation ---
    % Bernstein basis for degree 3
    b1 = 3 .* t .* (1-t).^2;
    b2 = 3 .* t.^2 .* (1-t);
    
    % --- Path Generation ---
    baseline = A' + t * dir_vec';
    deviation = (alpha1 .* b1 + alpha2 .* b2) * n_vec';
    P = baseline + deviation;
    
    % --- Cost Computation ---
    % Calling the separate function defined in the previous step
    [total_cost, details] = bernstein_cost(P, dt);
    
    % --- Visualization ---
    figure('Color', 'w');
    hold on; grid on; axis equal;
    
    % Plot the path
    plot(P(:,1), P(:,2), 'b-', 'LineWidth', 2.5, 'DisplayName', 'Drone Path');
    
    % Plot Start and End
    plot(A(1), A(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'Start (A)');
    plot(B(1), B(2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'End (B)');
    
    % Annotations
    xlabel('X Position'); ylabel('Y Position');
    title_str = sprintf('\\alpha_1 = %.2f, \\alpha_2 = %.2f', alpha1, alpha2);
    title(['Path Geometry: ', title_str]);
    
    % Cost Label (TextBox)
    dim = [.15 .6 .3 .3];
    str = {sprintf('Total Cost J: %.4f', total_cost), ...
           sprintf('Length Cost: %.4f', details.length), ...
           sprintf('Energy Cost: %.4f', details.energy)};
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
    
    legend('Location', 'northeast');
end

% --- LOCAL FUNCTION (As requested) ---
function [total_cost, details] = bernstein_cost(P, dt)
    v = diff(P) / dt;       % Velocity
    a = diff(v) / dt;       % Acceleration
    
    L_cost = sum(sqrt(sum(v.^2, 2))) * dt; % Arc Length
    E_cost = sum(sum(a.^2, 2)) * dt;       % Bending Energy
    
    total_cost = L_cost + 0.5 * E_cost; % Combined Cost
    details.length = L_cost;
    details.energy = E_cost;
end