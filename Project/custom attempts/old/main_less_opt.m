clc; clear; close all;
global Par

% --- Parameters ---
Par.A = [0; 0];      %start point
Par.B = [10; 0];     %end point
Par.dc = 0.01;       %time step for discretization
Par.LengthReference = norm(Par.B - Par.A);  % [m]

% Constraint parameters
Par.buffer           = 0.1;     % buffer value for obstacle constraints
Par.max_velocity     = 5.0;     % Maximum velocity constraint
Par.max_acceleration = 2.0;     % Maximum acceleration constraint
Par.max_curvature    = 0.5;     % Maximum curvature constraint


% --- obstacles --------------------------------------------------------------

obs_sparse     = [2,   0.2,  0.1;   
                  3,  -0.2,  0.1;     
                  4,    -1,  0.8;
                  6,     1,  0.8
                  9,     2,  1.5];  

num_per_side   = 10; 
gap            = 1.2; 
r              = 1;

%galley of obstacles
x_coords       = linspace(1, 9, num_per_side); 
obs_top        = [x_coords', ones(num_per_side, 1) * (gap/2 + r), ones(num_per_side, 1) * r];
obs_bottom     = [x_coords', ones(num_per_side, 1) * -(gap/2 + r), ones(num_per_side, 1) * r];
obs_gallery    = [obs_top; obs_bottom];

%diagonal line of obstacles
obs_pos        = [0 -1 r];
obs_diag       = [zeros(num_per_side,1), zeros(num_per_side,1),ones(num_per_side,1)*r];
obs_diag(1,:)  = obs_pos;
x_step         = 1;
y_step         = gap/(num_per_side);

for i=2:num_per_side
    obs_pos(1) = obs_pos(1) + x_step;
    obs_pos(2) = obs_pos(2) + y_step;
    obs_diag(i,:)  = obs_pos;
end

Par.obs = obs_gallery; % Combine all obstacles into one matrix


% --- Optimization Setup ---
n_vars   = 13;                                            % Change this to any integer!
x0       = [1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1]*2; 

cost     = cost_function(x0, Par)
[g, ceq] =            Constraint(x0, Par)