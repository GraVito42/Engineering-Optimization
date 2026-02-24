clc; clear; close all;
global Par

% --- Parameters ---
Par.A = [0; 0];
Par.B = [10; 0];
Par.dt = 0.01;

% --- (30 obstacles) ---
num_per_side = 10; 
x_coords = linspace(1, 9, num_per_side); 
gap = 1.2; r = 0.12;
obs_top = [x_coords', ones(num_per_side, 1) * (gap/2 + r), ones(num_per_side, 1) * r];
obs_bottom = [x_coords', ones(num_per_side, 1) * -(gap/2 + r), ones(num_per_side, 1) * r];
Par.obs = [obs_top ;obs_bottom];

%{
Par.obs = [2, 0.2, 0.1;   
           3, -0.2, 0.1;     
           4, -1, 0.8;
           6, 1, 0.8
           9, 2, 1.5];  
%}

% --- Optimization Setup ---
n_vars = 13; % Change this to any integer!
x0 = [1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1]*0.5; 
lb = ones(n_vars, 1) * -20;
ub = ones(n_vars, 1) * 20;

options = optimoptions('fmincon');
options.Display                     = 'iter-detailed';
options.Algorithm                   = 'sqp';
options.FunValCheck                 = 'off';       
options.MaxIter                     = 500;         
options.ScaleProblem                = true;        % Normalization of the variables
options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
options.FiniteDifferenceType        = 'central';
options.FiniteDifferenceStepSize    = 1e-2;
options.StepTolerance               = 1e-15; % Convergence criterion in step size
options.OptimalityTolerance         = 1e-9; % Convergence criterion in first order optimality
options.ConstraintTolerance         = 1e-4; % Determines the contraint tolerance
options.MaxFunEvals                 = 100000;
options.OutputFcn                   = {@TrajectoryPlotter};

[x_opt, FVAL] = fmincon(@(x) Objective(x, Par), x0, [], [], [], [], lb, ub, ...
    @(x) Constraint(x, Par), options);