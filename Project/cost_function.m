function total_cost = cost_function(alpha_n, Par)
% INPUTS:
%   P  : [Nx2] Matrix of trajectory points (x, y)
%   dx : Scalar, spatial step between samples
%
% OUTPUTS:
%   total_cost : Weighted sum of length and energy and safety
%   details    : Structure containing individual cost components

    global ub lb

    alpha = alpha_n .* (ub - lb) + lb;  % Scale from [0,1] to actual bounds

    % Compute curve features (Coords, Velocity, Acceleration, Curvature)
    P                   = bernstein_path(alpha,Par);
    [v, ~, j, k, dt, dx]    = kinematics(P, Par);


    % 1. Path length [m]  — sum of physical step lengths
    L1                   = sum(dx);
    L                   = L1 / Par.LengthReference - 1;  % Normalized length cost 
    
    % 2. Curvature Cost (Discrete Approximation)                           
    K1                   = sum(k.^2) * dt;             % Penalize squared curvature
    K                    = K1/Par.curvature_reference^2;     % Normalized curvature cost 

    % 3. Safety cost (Reciprocal of minimum distance to obstacles)
    [~, d]              = obstacle_distance(P, Par);
    Phi                 = abs(1/sum(d) - Par.buffer); 
    D                   = max(0, Phi);                      % Only penalize if within buffer distance

    % 4. Time cost
    v_n                 = vecnorm(v, 2,2);
    T1                   = sum(dx ./ v_n);
    T                    = T1 / (Par.LengthReference*Par.v_avg);  % Normalized time cost

    % 5. Jerk Cost
    j_mag_sq            = vecnorm(j, 2, 2).^2;                % (N-3) x 1
    j_ref_sq            = Par.jerk_reference^2;                       % align size with j
    J                   = sum(j_mag_sq ./ j_ref_sq) * dt;      % normalised jerk cost            

    total_weight = sum(Par.w);
    weights = Par.w / total_weight;                           % Normalize weights to sum to 1

    total_cost = weights(1)*L + weights(2)*K + weights(3)*D + weights(4)*T + weights(5)*J;

    % Return details for logging or plotting
    global details

    details.length      = L;
    details.curvature   = K;
    details.safety      = D;
    details.time        = T;
    details.jerk        = J;
end