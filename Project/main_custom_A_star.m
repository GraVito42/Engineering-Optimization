clc; clear; close all;
global Par ub lb var_history

addpath('island_images');

var_history = table();                                  % initialising debug arrays
hfig = figure('Name', 'Current Iteration results');     % initialising figure

% --- UAV model parameters ---
Par.A                = [800; -350];      %start point
Par.B                = [0; -400];        %end point
Par.v_avg            = 2.0;              %average velocity for our UAV
Par.dc               = 0.01;             % parametric step (t in [0,1]), controls curve resolution only

Par.LengthReference  = 1000;            % real physical distance A→B [m]

% --- Cost function weights ---
Par.w = [1.0 0.0 1.0 0.0 0.0];          % [Length Curvature Safety Time Jerk] cost weights

% Constraint parameters
Par.d_safe              = 5;            % Minimum safe distance from obstacle [m] // are we sure it's m?
Par.max_velocity        = 50;           % Maximum velocity constraint
Par.max_acceleration    = 20;           % Maximum acceleration constraint
Par.max_curvature       = 0.2;          % Maximum curvature constraint

% --- obstacles from the map image-----------------------------------------
[Par.obs, ~, ~, ~] = island_detection('../island_images/artic.png', false);

% --- using graph theory to find an initial path --------------------------
[raw_path, map] = find_global_path(Par);

% --- extract map bounds for variables ------------------------------------
min_x = map.LocalOriginInWorld(1);
min_y = map.LocalOriginInWorld(2);
max_x = min_x + map.GridSize(2);
max_y = min_y + map.GridSize(1);

% --- iterative optimisation setup ----------------------------------------
n_vars_start = 2;   
n_vars_max   = 30; 
converged     = false;
n_vars       = n_vars_start;

% --- Bernstein procedure initialisation ----------------------------------
dir_vec = Par.B - Par.A;
L = norm(dir_vec);
n_vec = [-dir_vec(2), dir_vec(1)] / L; % normal vector
vec_from_A = raw_path - Par.A(:)';         % distance between the found path and A
dist_along = (vec_from_A * dir_vec) / L; 

dc_raw = dist_along / L;               % found curve parametric step lenght

valid_idx = dc_raw > 0 & dc_raw < 1;   % filtering the non relevant points (L >1)
dc_raw = dc_raw(valid_idx);
vec_from_A = vec_from_A(valid_idx, :);

dev_target = vec_from_A * n_vec';

% --- setting fmincon Options ---------------------------------------------
options = optimoptions('fmincon');
options.Display                     = 'iter-detailed';
options.Algorithm                   = 'sqp';
options.FunValCheck                 = 'off';       
options.MaxIter                     = 500;         
options.ScaleProblem                = true;                     % Normalization of the variables
options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
options.FiniteDifferenceType        = 'central';
options.FiniteDifferenceStepSize    = 1e-3;                     % Why there was Par.d_safe?
options.StepTolerance               = 1e-15;                    % Convergence criterion in step size
options.OptimalityTolerance         = 1e-9;                     % Convergence criterion in first order optimality
options.ConstraintTolerance         = 1e-4;                     % Determines the contraint tolerance
options.MaxFunEvals                 = 100000;
options.OutputFcn                   = {@(x,optiomValues,state) TrajectoryPlotter_Astar(x,optiomValues,state,hfig, raw_path)};

tol_c = 5*options.ConstraintTolerance;
tol_fit = 12.0;                            % Max allowed geometric deviation from the found path

while n_vars <= n_vars_max && ~converged
    fprintf('\n--- Attempt with %d variables ---\n', n_vars);
    
    n_vars = 100;
    k = n_vars + 1;
    M = zeros(length(dc_raw), n_vars);   
    for i = 1:n_vars               % building bernstein polynomial
        coeff = nchoosek(k, i);
        M(:, i) = coeff .* (dc_raw.^i) .* ((1 - dc_raw).^(k - i));
    end

    % Least Squares resolution: Find optimal alphas to match the A* deviation
    x0_0      = M \ dev_target;
    dev_fit   = M*x0_0;
    fit_error = max(abs(dev_target - dev_fit));

    if fit_error > tol_fit
        fprintf('Level 1 FAILED at %d variables (Fit error: %.2f m > %.2f m).\n', n_vars, fit_error, tol_fit);
        fprintf('  Increasing n_vars.\n');
        n_vars = n_vars + 1; % Increase complexity
        continue; % Skip the rest and restart the loop
    else
        fprintf('  Level 1 PASSED at %d variables. Fit error: %.2f m.\n', n_vars, fit_error);
    end

    max_dev_allowed = max(max(abs(x0_0)) * 1.5, L/2);
    
    lb = -ones(n_vars, 1) * max_dev_allowed;
    ub =  ones(n_vars, 1) * max_dev_allowed;
    
    x0_norm = (x0_0 - lb) ./ (ub - lb);
    lb_n = zeros(n_vars, 1);
    ub_n = ones(n_vars, 1);

    [g, c_eq] = Constraint(x0_norm, Par);
    max_violation = max(g(:));
    if max_violation > tol_c
        fprintf('  Level 2 FAILED at %d variables (Max obstacle violation: %.4f > %.4f).\n', n_vars, max_violation, tol_c);
        fprintf('  Increasing n_vars.\n');
        n_vars = n_vars + 1; 
        continue; % Skip fmincon
    else
        fprintf('Level 2 PASSED at %d variables . Initial guess is safe.\n', n_vars);
    end
    
    [x_opt_norm, FVAL, exitflag, ~] = fmincon(@(x) cost_function(x, Par), x0_norm, ...
        [], [], [], [], lb_n, ub_n, ...
        @(x) Constraint(x, Par), options);
    
    if exitflag > 0
        fprintf('Optimization CONCLUDED successfully.\n');
        converged = true;
    else
        fprintf('Optimization failed. Increasing curve complexity...\n');
        n_vars = n_vars + 2; 
    end
end

if ~converged
    warning('fmincon failed to converge even with max variables (%d).', n_vars_max);
end

% Un-scale the final result
%X_opt_final = lb + x_opt_norm .* (ub - lb);


























% % --- Optimization Setup ---
% n_vars = 13;                                            % Change this to any integer!
% h      = norm(Par.B - Par.A) / (n_vars/2);  % Heuristic spacing for initial guess
% 
% %X0 = [1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1]*(h/2); 
% X0 = -ones(n_vars, 1) * h;  % Initial guess: straight line with uniform spacing
% lb = ones(n_vars, 1) *(-2)*h;
% ub = ones(n_vars, 1) * 2*h;
% 
% x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
% lb_n = zeros(n_vars,1);
% ub_n = ones(n_vars,1);
% 
% options = optimoptions('fmincon');
% options.Display                     = 'iter-detailed';
% options.Algorithm                   = 'sqp';
% options.FunValCheck                 = 'off';       
% options.MaxIter                     = 500;         
% options.ScaleProblem                = true;                     % Normalization of the variables
% options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
% options.FiniteDifferenceType        = 'central';
% options.FiniteDifferenceStepSize    = 1e-1;                     % Why there was Par.d_safe?
% options.StepTolerance               = 1e-15;                    % Convergence criterion in step size
% options.OptimalityTolerance         = 1e-9;                     % Convergence criterion in first order optimality
% options.ConstraintTolerance         = 1e-4;                     % Determines the contraint tolerance
% options.MaxFunEvals                 = 100000;
% options.OutputFcn                   = {@(x,optiomValues,state) TrajectoryPlotter(x,optiomValues,state,hfig)};
% 
% 
% [x_opt, FVAL] = fmincon(@(x) cost_function(x, Par), x0, [], [], [], [], lb, ub, ...
%     @(x) Constraint(x, Par), options);
