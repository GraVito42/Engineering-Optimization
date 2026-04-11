clc; clear; close all;
global Par ub lb

addpath('island_images');

% --- UAV model parameters ---
Par.A                = [0; 0];      %start point
Par.B                = [10; 0];     %end point
Par.v_avg            = 2.0;         %average velocity for our UAV
Par.dc               = 0.01;        % parametric step (t in [0,1]), controls curve resolution only

Par.LengthReference  = 100;  % real physical distance A→B [m]

% --- Cost function weights ---
Par.w = [1.0 1.0 10.0 1.0 1.0];     % [Length Curvature Safety Time Jerk] cost weights

% Constraint parameters
Par.d_safe              = 0.5;         % Minimum safe distance from obstacle [m]
Par.max_velocity        = 5.0;         % Maximum velocity constraint
Par.max_acceleration    = 3.0;         % Maximum acceleration constraint
Par.max_curvature       = 0.5;         % Maximum curvature constraint

% --- obstacles --------------------------------------------------------------

obs_sparse     = [2,   0.2,  0.1;   
                  3,  -0.2,  0.1;     
                  4,    -1,  0.8;
                  6,     1,  0.8
                  9,     2,  1.5];  

num_per_side   = 10; 
gap            = 1.2; 
r              = 0.3;

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

[Par.obs, ~, ~, ~] = island_detection('../island_images/canada.png', true);


% --- Optimization Setup ---
n_vars = 13;                                            % Change this to any integer!
h      = norm(Par.B - Par.A) / (n_vars/2);  % Heuristic spacing for initial guess

%X0 = [1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1]*(h/2); 
X0 = -ones(n_vars, 1) * h;  % Initial guess: straight line with uniform spacing
lb = ones(n_vars, 1) *(-2)*h;
ub = ones(n_vars, 1) * 2*h;

x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
lb_n = zeros(n_vars,1);
ub_n = ones(n_vars,1);

options = optimoptions('fmincon');
options.Display                     = 'iter-detailed';
options.Algorithm                   = 'sqp';
options.FunValCheck                 = 'off';       
options.MaxIter                     = 500;         
options.ScaleProblem                = true;                     % Normalization of the variables
options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
options.FiniteDifferenceType        = 'central';
options.FiniteDifferenceStepSize    = 1e-2;
options.StepTolerance               = 1e-15;                    % Convergence criterion in step size
options.OptimalityTolerance         = 1e-9;                     % Convergence criterion in first order optimality
options.ConstraintTolerance         = 1e-4;                     % Determines the contraint tolerance
options.MaxFunEvals                 = 100000;
options.OutputFcn                   = {@TrajectoryPlotter};


[x_opt, FVAL] = fmincon(@(x) cost_function(x, Par), x0, [], [], [], [], lb, ub, ...
    @(x) Constraint(x, Par), options);
