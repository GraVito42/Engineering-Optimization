clc; clear; close all;
global Par ub lb x0s idx 

addpath('island_images');

% --- UAV model parameters ---
Par.A                = [630; -400];      %start point
Par.B                = [0; -400];     %end point
Par.v_avg            = 2.0;         %average velocity for our UAV
Par.dc               = 0.01;        % parametric step (t in [0,1]), controls curve resolution only

Par.LengthReference  = 100;  % real physical distance A→B [m]

% --- Cost function weights ---
Par.w = [1.0 1.0 10.0 1.0 1.0];     % [Length Curvature Safety Time Jerk] cost weights

% Constraint parameters
Par.buffer              = 1/1000 * norm(Par.B - Par.A);         % buffer value for obstacle constraints
Par.max_velocity        = 5.0;         % Maximum velocity constraint
Par.max_acceleration    = 3.0;         % Maximum acceleration constraint
Par.max_curvature       = 0.5;         % Maximum curvature constraint
Par.curvature_reference = 0.1; 
Par.jerk_reference      = 1;          % Reference value to normalise jerk

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

[Par.obs, ~, ~, ~] = island_detection('../island_images/artic.png', false);


% --- Optimization Setup ---
n_vars = 13;                                            % Change this to any integer!
h      = norm(Par.B - Par.A) / (n_vars/2);  % Heuristic spacing for initial guess

%X0 = [1, 1, 1, 1, 1, 1, 0, -1, -1, -1, -1, -1, -1]*(h/2); 
lb = ones(n_vars, 1) *(-2)*h;
ub = ones(n_vars, 1) * 2*h;

slices = 10;
fc  = zeros(4*slices+2,3);   
x0s = zeros(4*slices+2, n_vars);


for j=-slices:slices
    H = j*(2*h/slices);
    n_fc = j+slices+1;

    X0 = ones(n_vars, 1) * H;  % Initial guess: straight line with uniform spacing
    x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
    fc(n_fc,:) = [cost_function(x0, Par)/max(Constraint(x0, Par)), max(Constraint(x0, Par)), n_fc];
    x0s(n_fc,:) = x0;
end

for j=-slices:slices
    H = j*(2*h/slices);
    n_fc = j+3*slices+2;

    X01 = ones((n_vars-1)/2, 1);
    X02 = -ones((n_vars-1)/2, 1);
    X0 = [X01; 0; X02] * H;  %
    x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
    fc(n_fc,:) = [cost_function(x0, Par)/max(Constraint(x0, Par)), max(Constraint(x0, Par)), n_fc];
    x0s(n_fc,:) = x0;
end

fcs = fc(fc(:, 2) < 1e-4, :);
[fcmin, idxs] = min(fcs(:, 1));


if isempty(fcmin)
    
    % First systematic trial
    fc_extra = zeros(n_vars*slices,3);
    x0s_extra = zeros(n_vars*slices, n_vars);

    % Search for S paths with a moving focal point
    for k=1:n_vars

            % First case with the focal on the first variable
            if k == 1
                for j=-slices:slices
                    H = j*(2*h/slices);
                    n_fc = j+slices+1;
    
                    X0 = [0; ones(n_vars-1, 1)] * H;  %
                    x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
                    fc_extra(n_fc,:) = [cost_function(x0, Par)/max(Constraint(x0, Par)), max(Constraint(x0, Par)), n_fc];
                    x0s_extra(n_fc,:) = x0;
                end
            % Last case with the focal point on the last variable
            elseif k == n_vars
                for j=-slices:slices
                    H = j*(2*h/slices);
                    n_fc = j+(2*k-1)*slices+k;
    
                    X0 = [ones(n_vars-1, 1); 0] * H;  %
                    x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
                    fc_extra(n_fc,:) = [cost_function(x0, Par)/max(Constraint(x0, Par)), max(Constraint(x0, Par)), n_fc];
                    x0s_extra(n_fc,:) = x0;
                end
            % Generic case with the focal point in position 2 to n_vars - 1
            else
                for j=-slices:slices
                    H = j*(2*h/slices);
                    n_fc = j+(2*k-1)*slices+k;
    
                    X01 = ones(k-1, 1);
                    X02 = -ones(n_vars-k, 1);
                    X0 = [X01; 0; X02] * H;  %
                    x0 = (X0 - lb)./(ub - lb);  % Normalise initial guess to [0,1]
                    fc_extra(n_fc,:) = [cost_function(x0, Par)/max(Constraint(x0, Par)), max(Constraint(x0, Par)), n_fc];
                    x0s_extra(n_fc,:) = x0;
                end
            end
          
    end
    
    % Update the checked paths arrays
    x0s = [x0s; x0s_extra];
    fc  = [fc; fc_extra];

    % Extract from the evaluated function the feasible ones, and pick the
    % smallest f_val/constraint
    fcs_extra = fc_extra(fc_extra(:, 2) < 1e-4, :);
    [fcmin_extra, idx_extra] = min(fcs_extra(:, 1));
    n   = length(x0s);
    
    if ~isempty(fcmin_extra)
        idx = 4*slices+2 + idx_extra;
        
    else

        %Second trial: random sampling
        check     = 1;
        step      = 0;
        step_step = 0;
        n_steps   = 10000;
        fc_rand  = ones(check, 3);
    
        %Random search for the feasible starting point
        while (check < 6 && step < n_steps) 
    
            step = step + 1;
            x0 = rand(n_vars,1);
            [c, tilde] = Constraint(x0, Par);
    
            % If a point is feasible, add it to the array
            if max(c) < 1e-4
                fc_rand(check,:) = [cost_function(x0, Par)/max(c), max(c), check];
                x0s(n+check,:) = x0;
                check = check + 1;
                fprintf('Found %d/%d feasible initial guesses after %d random samples.\n', ...
                        length(fc_rand)-check, length(fc_rand), step);
            else
                if mod(step, 1000) == 0
                fprintf('Checked %d random samples, found %d feasible initial guesses so far.\n', ...
                        step, check-1);
                end
            end
    
            %Check every n_steps to see if we found something
            if (step == n_steps && check == 1)
                step = 0;
                fprintf('No random feasible point found, looping again .\n');
                step_step = step_step +1;
                if step_step == 10
                    step = n_steps;
                end
            end
            
        end

        %Pick the feasible point found with the lowest function value
        [fcmin_rand, idx_rand] = min(fc_rand(:, 1));
        idx = n + idx_rand;



        if check == 1

            %Third way: if nothing else functioned,
            % take the smallest constraint violation
            [fcmin, idx] = min(fc(:, 2));
        end

     end

else
    idx = fcs(idxs,3);
end

lb_n = zeros(n_vars, 1);   % normalized lower bound
ub_n = ones(n_vars, 1);    % normalized upper bound
%x0 = x0s(idx, :);
%
%Plotinitialguess(Par, x0s, idx);
% 
% options = optimoptions('fmincon');
% options.Display                     = 'iter-detailed';
% options.Algorithm                   = 'sqp';
% options.FunValCheck                 = 'off';       
% options.MaxIter                     = 500;         
% options.ScaleProblem                = true;                     % Normalization of the variables
% options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
% options.FiniteDifferenceType        = 'central';
% options.FiniteDifferenceStepSize    = 1e-3;
% options.StepTolerance               = 1e-15;                    % Convergence criterion in step size
% options.OptimalityTolerance         = 1e-9;                     % Convergence criterion in first order optimality
% options.ConstraintTolerance         = 1e-4;                     % Determines the contraint tolerance
% options.MaxFunEvals                 = 100000;
% options.OutputFcn                   = {@TrajectoryPlotter_custom};
% 
% [x_opt, FVAL] = fmincon(@(x) cost_function(x, Par), x0, [], [], [], [], lb, ub, ...
%     @(x) Constraint(x, Par), options);

x_opts = zeros(idx_rand, n_vars);
fvals  = zeros(idx_rand, 1);
c_opts = zeros(idx_rand, 1);

for i=1:check-1
    x0 = x0s(n+i,:);

    options = optimoptions('fmincon');
    options.Display                     = 'iter-detailed';
    options.Algorithm                   = 'sqp';
    options.FunValCheck                 = 'off';       
    options.MaxIter                     = 500;         
    options.ScaleProblem                = true;                     % Normalization of the variables
    options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
    options.FiniteDifferenceType        = 'central';
    options.FiniteDifferenceStepSize    = 1e-3;
    options.StepTolerance               = 1e-15;                    % Convergence criterion in step size
    options.OptimalityTolerance         = 1e-9;                     % Convergence criterion in first order optimality
    options.ConstraintTolerance         = 1e-4;                     % Determines the contraint tolerance
    options.MaxFunEvals                 = 100000;
    options.OutputFcn                   = {@TrajectoryPlotter_custom};

    [x_opt, FVAL] = fmincon(@(x) cost_function(x, Par), x0, [], [], [], [], lb, ub, ...
        @(x) Constraint(x, Par), options);

    x_opts(i,:) = x_opt;
    fvals(i)    = FVAL;
    c_opts(i)     = max(Constraint(x_opt, Par));
end

% if min(c_opts) > 1e-4
%     every_step = 0:0.5:1;
%     every      = cell(1, n_vars);
%     [every{:}] = ndgrid(every_step);
%     every_comb = cell2mat(cellfun(@(x) x(:), every, 'UniformOutput', false));
    
%     x_opts_ev = zeros(numel(every_comb), n_vars);
%     fvals_ev  = zeros(numel(every_comb), 1);
%     c_opts_ev = zeros(numel(every_comb), 1);
%     counter   = 0;

%     for i=1:numel(every_comb)
%         counter = counter +1;
%         x0 = every_comb(i);
        
%         c_opts_ev(i) = max(Constraint(x0,Par));

%         if mod(counter, 10000) == 0
%                 fprintf('Checked %d samples \n', ...
%                         counter);
%         end

%         if c_opts_ev(i) < 0
%             fprintf("FOUND FEASIBLE!")

%             options = optimoptions('fmincon');
%             options.Display                     = 'iter-detailed';
%             options.Algorithm                   = 'sqp';
%             options.FunValCheck                 = 'off';       
%             options.MaxIter                     = 500;         
%             options.ScaleProblem                = true;                     % Normalization of the variables
%             options.PlotFcn                     = {@optimplotfval, @optimplotx, @optimplotfirstorderopt, @optimplotstepsize, @optimplotconstrviolation, @optimplotfunccount};
%             options.FiniteDifferenceType        = 'central';
%             options.FiniteDifferenceStepSize    = 1e-3;
%             options.StepTolerance               = 1e-15;                    % Convergence criterion in step size
%             options.OptimalityTolerance         = 1e-9;                     % Convergence criterion in first order optimality
%             options.ConstraintTolerance         = 1e-4;                     % Determines the contraint tolerance
%             options.MaxFunEvals                 = 100000;
%             options.OutputFcn                   = {@TrajectoryPlotter_custom};
        
%             [x_opt, FVAL] = fmincon(@(x) cost_function(x, Par), x0, [], [], [], [], lb, ub, ...
%                 @(x) Constraint(x, Par), options);
        
%             x_opts(i,:) = x_opt;
%             fvals(i)    = FVAL;

%             break
%         end
%     end
% end

