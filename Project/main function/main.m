clc; clear; close all;
global Par ub lb var_history

addpath('island_images');

var_history = table();                                  % initialising debug arrays

% --- UAV model parameters ---
Par.A                = [0;       100];      % start point
Par.B                = [1100;   -450];      % end point
Par.v_avg            = 10;              %average velocity for our UAV
Par.dc               = 1e-3;             % parametric step (t in [0,1]), controls curve resolution only

%Par.LengthReference  = 1000;            % real physical distance A→B [m]
Par.LengthReference = norm(Par.B - Par.A);

% --- Cost function weights ---
Par.w = [1.0 1.0 0.2 0 1.0];          % [Length Curvature Safety Time Jerk] cost weights

% Constraint parameters
Par.d_safe              = 5;            % Minimum safe distance from obstacle [-]
Par.max_velocity        = 50;           % Maximum velocity constraint
Par.max_acceleration    = 20;           % Maximum acceleration constraint
Par.max_curvature       = 0.2;          % Maximum curvature constraint

% --- obstacles --------------------------------------------------------------
[Par.obs, ~, ~, ~] = island_detection('../island_images/canada.png', false);


% --- Optimization Setup ---
bernstein_coeff = 13;                                     % Bernstein order -1
n_vars          = 2 * bernstein_coeff;                                            
bound_coeff     = 3;                                      % Defines bounds dimension
h               = norm(Par.B - Par.A) / bernstein_coeff;  % Heuristic spacing for initial guess


%X0 = [1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1]*(h/2); 
X0 = zeros(n_vars, 1);  % Initial guess: straight line with uniform spacing

% Bounds
% lb = ones(n_vars, 1) *(-bound_coeff)*h;    % Old used bounds, here for backup
% ub = ones(n_vars, 1) * bound_coeff*h;

% Normal: symmetric, free to deviate left or right
lb_eta   = -bound_coeff * h * ones(bernstein_coeff,1);
ub_eta   =  bound_coeff * h * ones(bernstein_coeff,1);

% Tangential: asymmetric — allow compression but forbid fold-back
% The critical limit is -h (one node spacing), stay within 90% of it
lb_tau   = -0.9 * h * ones(bernstein_coeff,1);
ub_tau   =  bound_coeff * h * ones(bernstein_coeff,1);

lb = [lb_eta; lb_tau];
ub = [ub_eta; ub_tau];

x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
lb_n = zeros(n_vars,1);
ub_n = ones(n_vars,1);

hfig = figure('Name', 'Current Iteration results');             % initialising figure
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
options.OutputFcn                   = {@(x,optiomValues,state) TrajectoryPlotter(x,optiomValues,state,hfig)};

 % num_starts = 2;
 % [x_opt, FVAL] = multi_start_optimization(Par, n_vars, num_starts, options, lb_n, ub_n);

[x_opt, FVAL] = fmincon(@(x) cost_function(x, Par), x0, [], [], [], [], lb_n, ub_n, ...
  @(x) Constraint(x, Par), options);
