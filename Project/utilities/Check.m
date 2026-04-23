%% Post-Optimization Consistency Check
% Run this script after fmincon finishes to verify the physical results.

addpath('..\main function\')

fprintf('--- Optimization Consistency Report ---\n');

% 1. Reconstruct physical trajectory and kinematics
global ub lb 
alpha_opt = x_opt .* (ub - lb) + lb;
P_opt = bernstein_path(alpha_opt, Par);
[v, a, j, k, dt, dx] = kinematics(P_opt, Par);

% Magnitudes
v_mag = vecnorm(v, 2, 2);
a_mag = vecnorm(a, 2, 2);
j_mag = vecnorm(j, 2, 2);
L_total = sum(dx);
T_total = sum(dx ./ v_mag);

% 2. Check Constraint Saturation (Are we within limits?)
max_v_calc = max(v_mag);
max_a_calc = max(a_mag);
max_k_calc = max(k);

fprintf('\n1. KINEMATIC LIMITS:\n');
fprintf('   Velocity:     Max = %.2f m/s (Limit: %.2f) | Status: %s\n', ...
    max_v_calc, Par.max_velocity, check_limit(max_v_calc, Par.max_velocity));
fprintf('   Acceleration: Max = %.2f m/s² (Limit: %.2f) | Status: %s\n', ...
    max_a_calc, Par.max_acceleration, check_limit(max_a_calc, Par.max_acceleration));
fprintf('   Curvature:    Max = %.4f 1/m (Limit: %.2f) | Status: %s\n', ...
    max_k_calc, Par.max_curvature, check_limit(max_k_calc, Par.max_curvature));


%Obstacle Safety Verification
[~, d_vec] = obstacle_distance(P_opt, Par);
min_d = min(d_vec);
collision_points = sum(d_vec <= 0);
buffer_violations = sum(d_vec > 0 & d_vec < Par.d_safe);

fprintf('\n2. OBSTACLE SAFETY:\n');
fprintf('   Minimum Clearance:  %.2f m (Safe Buffer: %.2f m)\n', min_d, Par.d_safe);
if min_d <= 0
    fprintf('   STATUS:             FAILED (Collision detected at %d points)\n', collision_points);
elseif min_d < Par.d_safe
    fprintf('   STATUS:             CAUTION (Inside buffer at %d points)\n', buffer_violations);
else
    fprintf('   STATUS:             CLEAR\n');
end


% Cost Normalization Breakdown (Weights vs. Values)
cost_function(x_opt, Par);
global details
w_norm = Par.w / sum(Par.w); % Normalized weights used in cost_function

fprintf('\n3. COST ANALYSIS (Weighting Effectiveness):\n');
fprintf('   Component   | Norm. Value | Weight | Weighted Contrib.\n');
fprintf('   ------------|-------------|--------|------------------\n');
fprintf('   Length (L)  |   %8.4f  |  %.2f  |   %.4f\n', details.length, w_norm(1), details.length * w_norm(1));
fprintf('   Curvature(K)|   %8.4f  |  %.2f  |   %.4f\n', details.curvature, w_norm(2), details.curvature * w_norm(2));
fprintf('   Safety (D)  |   %8.4f  |  %.2f  |   %.4f\n', details.safety, w_norm(3), details.safety * w_norm(3));
fprintf('   Time (T)    |   %8.4f  |  %.2f  |   %.4f\n', details.time, w_norm(4), details.time * w_norm(4));
fprintf('   Jerk (J)    |   %8.4f  |  %.2f  |   %.4f\n', details.jerk, w_norm(5), details.jerk * w_norm(5));


% Geometric & Time Sanity
T_ref = Par.LengthReference / Par.v_avg;
dist_AB = norm(Par.B - Par.A);

fprintf('\n3. PHYSICAL SANITY:\n');
fprintf('   Total Distance: %.2f m (Reference A->B: %.2f m)\n', L_total, dist_AB);
fprintf('   Flight Time:    %.2f s (Ref Time: %.2f s)\n', T_total, T_ref);
fprintf('   Avg Velocity:   %.2f m/s (Target: %.2f m/s)\n', L_total/T_total, Par.v_avg);


% Physical Sanity & Smoothness
T_ref = Par.LengthReference / Par.v_avg;
fprintf('\n4. PHYSICAL SANITY & SMOOTHNESS:\n');
fprintf('   Flight Time:    %.2f s (Ref: %.2f s) | Variance: %.1f%%\n', T_total, T_ref, (T_total/T_ref-1)*100);
fprintf('   Avg Jerk:       %.2f m/s³\n', mean(j_mag));
fprintf('   Max Jerk:       %.2f m/s³\n', max(j_mag));


% Variable Bound Check (Is the search space too small?)
at_ub = sum(x_opt > 0.99);
at_lb = sum(x_opt < 0.01);
fprintf('\n4. SEARCH SPACE:\n');
fprintf('   Variables at Upper Bound: %d / %d\n', at_ub, length(x_opt));
fprintf('   Variables at Lower Bound: %d / %d\n', at_lb, length(x_opt));

if at_ub + at_lb > 0
    warning('Optimization hit search bounds. Path might be "clipped". Consider increasing lb/ub.');
end

function str = check_limit(val, limit)
    if val <= limit * 1.001 % 0.1% tolerance
        str = 'OK';
    else
        str = 'VIOLATED';
    end
end