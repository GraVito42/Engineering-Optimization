clc; clear; close all;

global Par
% --- Parameters ---
Par.A = [0; 0];
Par.B = [10; 0];
Par.dt = 0.01;

% Obstacles: [x, y, radius]

%{
num_per_side = 10; 
x_coords = linspace(1, 9, num_per_side); 
gap = 1.2;
r = 0.34;   

obs_top = [x_coords', ones(num_per_side, 1) * (gap/2 + r), ones(num_per_side, 1) * r];
obs_bottom = [x_coords', ones(num_per_side, 1) * -(gap/2 + r), ones(num_per_side, 1) * r];

Par.obs = [obs_top; obs_bottom];
%}

obs = zeros(5,3);

for i=1:5
    k = 10 * rand();
    j = -5 + (5 +5) * rand();
    r = 0.1 + (1 - 0.1) * rand();
    obs(i,:) = [k j r];
end

Par.obs = obs;

%{
Par.obs = [2, 0.2, 0.1;   
           3, -0.2, 0.1;     
           4, -1, 0.8;
           6, 1, 0.8
           9, 2, 1.5];  
%}

% --- Optimization Setup ---
lb = [-20, -20];
ub = [20, 20];
x0 = [2, -2]; 

options = optimoptions('fmincon');
options.Display                     = 'iter-detailed';
options.Algorithm                   = 'sqp';
options.FunValCheck                 = 'off';       
options.MaxIter                     = 500;         
options.ScaleProblem                = true;        % Normalization of the variables
options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
options.FiniteDifferenceType        = 'central';
options.FiniteDifferenceStepSize    = 1e-2;
options.StepTolerance               = 1e-9; % Convergence criterion in step size
options.OptimalityTolerance         = 1e-9; % Convergence criterion in first order optimality
options.ConstraintTolerance         = 1e-4; % Determines the contraint tolerance
options.MaxFunEvals                 = 100000;
options.OutputFcn                   = {@TrajectoryPlotter};

tic;
[x_opt, FVAL, EXITFLAG] = fmincon(@(x) Objective(x, Par), x0, [], [], [], [], lb, ub, ...
    @(x) Constraint(x, Par), options);
toc;

fprintf('Optimized Alphas: [%.2f, %.2f]\n', x_opt(1), x_opt(2));