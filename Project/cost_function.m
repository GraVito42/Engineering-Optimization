function total_cost = cost_function(alpha, Par)
% INPUTS:
%   P  : [Nx2] Matrix of trajectory points (x, y)
%   dx : Scalar, spatial step between samples
%
% OUTPUTS:
%   total_cost : Weighted sum of length, energy, safety, time and jerk
%   details    : Structure containing individual cost components
%
% NORMALIZATION STRATEGY (hybrid):
%   L, T  — normalized by straight-line reference (always well-defined, ≈1 for direct path)
%   K     — normalized by max_curvature^2 * T_ref: = 1 if path is at curvature limit everywhere
%   J     — normalized by j_max^2 * T_ref where j_max = max_acceleration * v_avg / LengthReference
%   D     — physical barrier formulation with d_safe [m] (geometric meaning preserved)
    
    % Compute curve features (Coords, Velocity, Acceleration, Curvature)
    P                       = bernstein_path(alpha,Par);
    [v, ~, j, k, dt, dx]    = kinematics(P, Par);
    v_n                     = vecnorm(v, 2, 2);
    Tref                    = Par.LengthReference/Par.v_avg; % Time to fly straight from A to B [s]

    % 1. Length cost
    % Normalized by straight-line distance A→B.
    % c_L = 1 for a straight path, > 1 for any deviation.
    L1                  = sum(dx); % True arc length
    L                   = L1 / Par.LengthReference;  % Normalized length cost ≈ 1 for near-straight path
    
    % 2. Curvature cost
    % Integral of squared curvature over time (bending energy).
    % Normalized by value at initial guess so c_K = 1 at x0.
    K1  = sum(k.^2) * dt;                       % [1/m² · s]
    K   = K1 / (Par.max_curvature^2 * Tref);   % = 1 at x0 [.]

    % 3. Safety cost (Reciprocal of minimum distance to obstacles)
    % Barrier penalty: zero when d_min >= d_safe, grows as path approaches obstacles.
    % d_safe has clear physical meaning [m]; no dependence on obstacle count.
    [~, d]              = obstacle_distance(P, Par);
    phi                 = max(0, 1/min(d) - 1/Par.d_safe);  % Only penalize if within buffer distance
    D                   = (phi * Par.d_safe)^2;             % Dimensionless and = 0 when safe

    % 4. Time cost
    T1                   = sum(dx ./ v_n);  % actual travel time [s]
    T                    = T1 / Tref;       % Normalized time cost ≈ 1 for near-straight path [.]

    % 5. Jerk Cost
    % Integral of squared jerk magnitude over time (smoothness).
    % Normalized by value at initial guess so c_J = 1 at x0.
    J1  = sum(vecnorm(j, 2, 2).^2) * dt;                                % [m²/s⁶ · s]
    j_max   = Par.max_acceleration * Par.v_avg / Par.LengthReference;   % max jerk from kinematic limits [m/s³]
    J   = J1 / (j_max^2 * Tref);                                       % = 1 at x0 [.]          

    % Weighted sum
    weights = Par.w / sum(Par.w);  % Normalize weights to sum to 1
    total_cost = weights(1)*L + weights(2)*K + weights(3)*D + weights(4)*T + weights(5)*J;

    % Return details for logging or plotting
    global details

    details.length      = L;
    details.curvature   = K;
    details.safety      = D;
    details.time        = T;
    details.jerk        = J;
end