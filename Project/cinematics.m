function [v, a, j, k, dt] = cinematics(P, Par)

    L            = sum(vecnorm(diff(P), 2, 2));    % Total path length
    T            = L / Par.v_avg;                  % Total time based on average velocity
    dt           = T*Par.dx;                       % Spatial step for finite differences

    % 1. Compute Derivatives (Velocity, Acceleration, Curvature)
    v = diff(P) / dt;               % Velocity vectors m   [ (N-1) x 2 ]
    a = diff(v) / dt;               % Acceleration vectors [ (N-2) x 2 ]
    j = diff(a) / dt;               % Jerk vectors         [           ]
    
    % Curvature k = |v x a| / |v|^3
    v_mid = v(1:end-1,:) + v(2:end,:); % Midpoint velocity for curvature
    k = vecnorm(v_mid(:,1).*a(:,2) - v_mid(:,2).*a(:,1), 2, 2) ...
        ./ (vecnorm(v_mid, 2, 2).^3 + 1e-6); 

end