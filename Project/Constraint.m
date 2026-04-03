function [g, ceq] = Constraint(alpha, Par)

    P = bernstein_path(alpha,Par);

    [g, ~] = obstacle_distance(P, Par); % Obstacle constraint

    [v, a, ~, k] = cinematics(P, Par);

    max_velocity = Par.max_velocity;
    v_mag = vecnorm(v, 2, 2); 
    g = [g; v_mag - max_velocity];      % Velocity constraint
    
    max_acceleration = Par.max_acceleration;
    a_mag = vecnorm(a, 2, 2);  
    g = [g; a_mag - max_acceleration];  % Acceleration constraint
    
    max_curvature = Par.max_curvature;
    g = [g; k - max_curvature];         % Curvature constraint

    ceq = [];
end