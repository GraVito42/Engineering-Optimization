%% Post-Optimization Consistency Check
% Run this script after fmincon finishes to verify the physical results.

fprintf('--- Optimization Consistency Report ---\n');

% 1. Reconstruct physical trajectory and kinematics
global ub lb 
alpha_opt = x_opt .* (ub - lb) + lb;
P_opt = bernstein_path(alpha_opt, Par);
[v, a, j, k, dt, dx] = kinematics(P_opt, Par);

% Magnitudes
v_mag = vecnorm(v, 2, 2);
a_mag = vecnorm(a, 2, 2);
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

% 3. Check Normalization Scaling (Are cost components balanced?)
% Run cost function once more to populate global 'details'
cost_function(x_opt, Par);
global details

fprintf('\n2. COST NORMALIZATION (Should be O(1)):\n');
fprintf('   Length (L):    %.4f\n', details.length);
fprintf('   Curvature (K): %.4f\n', details.curvature);
fprintf('   Safety (D):    %.4f\n', details.safety);
fprintf('   Time (T):      %.4f\n', details.time);
fprintf('   Jerk (J):      %.4f\n', details.jerk);

% 4. Geometric & Time Sanity
T_ref = Par.LengthReference / Par.v_avg;
dist_AB = norm(Par.B - Par.A);

fprintf('\n3. PHYSICAL SANITY:\n');
fprintf('   Total Distance: %.2f m (Reference A->B: %.2f m)\n', L_total, dist_AB);
fprintf('   Flight Time:    %.2f s (Ref Time: %.2f s)\n', T_total, T_ref);
fprintf('   Avg Velocity:   %.2f m/s (Target: %.2f m/s)\n', L_total/T_total, Par.v_avg);

% 5. Variable Bound Check (Is the search space too small?)
at_ub = sum(x_opt > 0.99);
at_lb = sum(x_opt < 0.01);
fprintf('\n4. SEARCH SPACE:\n');
fprintf('   Variables at Upper Bound: %d / %d\n', at_ub, length(x_opt));
fprintf('   Variables at Lower Bound: %d / %d\n', at_lb, length(x_opt));
if at_ub + at_lb > 0
    warning('Optimizer hit bounds. Consider increasing lb/ub range.');
end

function str = check_limit(val, limit)
    if val <= limit * 1.001 % 0.1% tolerance
        str = 'OK';
    else
        str = 'VIOLATED';
    end
end