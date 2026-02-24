function [cost] = bernstein_path(alpha1, alpha2)
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
    cost = details.length;
end
