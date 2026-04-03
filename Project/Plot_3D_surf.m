clc
clear all
close all

% --- Parameters ---
Par.A = [0; 0];
Par.B = [10; 10];
Par.dx = 0.1;

% --- Grid Setup ---
res = 100;
l = 20;
alpha1_vec = linspace(-l, l, res);
alpha2_vec = linspace(-l, l, res);

% Pre-allocate Z
Z = zeros(res, res);

% --- Computation ---
for i = 1:res     % Iterate over alpha2 (Y-axis/Rows)
    for j = 1:res % Iterate over alpha1 (X-axis/Cols)
        % Note the order: j corresponds to alpha1, i to alpha2
        alpha = [alpha1_vec(j), alpha2_vec(i)]; 
        
        P = bernstein_path(alpha, Par);
        
        % Using the correct function name from your files
        [cost, ~] = cost_function(P, Par.dx); 
        
        % Assign to Z(row, col) -> Z(y_index, x_index)
        Z(i,j) = cost;
    end
end

% --- Visualization ---
figure('Color', 'w');
surf(alpha1_vec, alpha2_vec, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.8);

% Aesthetics
colormap jet;
colorbar;
xlabel('\alpha_1');
ylabel('\alpha_2');
zlabel('Cost');
title('Cost Landscape (Squared Curvature + Length)');
view(-45, 30);
grid on;