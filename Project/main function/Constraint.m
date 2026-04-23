function [g, ceq] = Constraint(alpha_n, Par)

    global ub lb

    alpha = alpha_n .* (ub - lb) + lb;  % Scale from [0,1] to actual bounds
    P = bernstein_path(alpha, Par);
    [v, a, ~, k, ~, ~] = kinematics(P, Par);

    % 1. Obstacle constraint
    [g_obs, ~] = obstacle_distance(P, Par); % (r + d_safe) - min_dist [m]
    g = g_obs ./ Par.LengthReference;

    % 2. Velocity constraint
    v_mag = vecnorm(v, 2, 2);
    g = [g; (v_mag - Par.max_velocity) / Par.max_velocity];

    % 3. Acceleration constraint
    a_mag = vecnorm(a, 2, 2);
    g = [g; (a_mag - Par.max_acceleration) / Par.max_acceleration];

    % 4. Curvature constraint
    g = [g; (k - Par.max_curvature) / Par.max_curvature];

    ceq = [];
end