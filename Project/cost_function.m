function total_cost = cost_function(alpha_n, Par)
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
    global lb ub
    alpha = alpha_n .* (ub - lb) + lb;
    P                       = bernstein_path(alpha, Par);
    [v, ~, j, k, dt, dx]    = kinematics(P, Par);
    v_n                     = vecnorm(v, 2, 2);
    Tref                    = Par.LengthReference/Par.v_avg; % Time to fly straight from A to B [s]

    % kinematic references 
    P_star   = bernstein_path(ub, Par);
    [~, ~, j_ref, k_ref, dt_ref, ~]    = kinematics(P_star, Par);

    % 1. Length cost
    % Normalized by straight-line distance A→B.
    % c_L = 1 for a straight path, > 1 for any deviation.
    L1                  = sum(dx); % True arc length
    L                   = L1 / Par.LengthReference;  % Normalized length cost ≈ 1 for near-straight path
    
    % 2. Curvature cost
    % Physical normalization: Integral of (max_curvature^2) over Tref
    K1  = sum(k.^2) * dt;                       %           [1/m² · s]
    K_limit_integral = (Par.max_curvature^2) * Tref;
    %K   = K1 / (sum(k_ref.^2) * dt_ref);    % = 1 at x0 [.] ????
    K   = K1 / K_limit_integral;

    % 3. Safety cost (Reciprocal of minimum distance to obstacles)
    % Barrier penalty: zero when d_min >= d_safe, grows as path approaches obstacles.
    % d_safe has clear physical meaning [m]; no dependence on obstacle count.
    [~, d]              = obstacle_distance(P, Par);
    if min(d) <= 0
        % path is inside obstacle, so apply a penalty proportional to depth
        %D    = 1 - min(d)^2/(Par.LengthReference * Par.d_safe); % >1
        D    = (1 - min(d)/(Par.d_safe)); % >1 
    else
        D    = max(0, 1 - min(d)/(Par.d_safe));  % Only penalize if within buffer distance
    end

    % 4. Time cost
    T1                   = sum(dx ./ v_n);  % actual travel time [s]
    T                    = T1 / Tref;       % Normalized time cost ≈ 1 for near-straight path [.]

    % 5. Jerk Cost
    % Define a physical jerk scale: how fast we can realistically change acceleration
    % Heuristic: max_acceleration / (time to reach avg velocity)
     J1  = sum(vecnorm(j, 2, 2).^2) * dt;                                % [m²/s⁶ · s]
     j_scale = Par.max_acceleration / Tref;
     J_limit_integral = (j_scale^2) * Tref;
     J   = J1 / J_limit_integral;

    
    % % ----- Define Physical Jerk Scale ----- doesn't seem to change much
    % % from the one above honestly 
    % % We define a "jerk limit" as the ability to reach max acceleration 
    % % in the time it takes to travel the reference distance.
    % % j_limit = Par.max_acceleration / Tref; 
    % % 
    % % 3. Calculate Normalized Jerk Cost
    % % Integral of ||j||^2 dt
    % % j_mag_sq = sum(j.^2, 2); % Squared magnitude of jerk vector
    % % 
    % % We must match the length of dt_vec to the length of j (N-3)
    % % Use the mid-point time steps calculated in kinematics
    % % dt_accel = 0.5 * (dt_vec(1:end-1) + dt_vec(2:end));
    % % dt_jerk_vec = 0.5 * (dt_accel(1:end-1) + dt_accel(2:end));
    % % 
    % % Numerical Integral: Sum(j^2 * dt)
    % % J_total = sum(j_mag_sq) * dt; % scalar .* dt_jerk_vec);
    % % 
    % % Physical Normalization Factor (Max possible jerk integrated over Tref)
    % % J_norm_factor = (j_limit^2) * Tref;
    % % 
    % % Final Normalized Jerk Cost
    % % J = J_total / J_norm_factor;


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