function total_cost = cost_function(alpha, Par)
% INPUTS:
%   P  : [Nx2] Matrix of trajectory points (x, y)
%   dx : Scalar, spatial step between samples
%
% OUTPUTS:
%   total_cost : Weighted sum of length and energy and safety
%   details    : Structure containing individual cost components

    % Compute curve features (Coords, Velocity, Acceleration, Curvature)
    P                   = bernstein_path(alpha,Par);
    [v, ~, j, k, dt]    = cinematics(P, Par);

    % 1. Distance Cost (L2 Norm of velocity)
    % Sum of infinitesimal displacements: ds = |v| * dt
    step_lengths        = norm(v);
    L                   = sum(step_lengths);
    
    % 2. Curvature Cost (Discrete Approximation)                           
    K                   = norm(k)^2 * dt;             % Penalize squared curvature
    
    % 3. Safety cost (Reciprocal of minimum distance to obstacles)
    [~, d]              = obstacle_distance(P, Par);
    Phi                 = (1/sum(d) - 1/Par.buffer)^2; 
    D                   = max(0, Phi);                      % Only penalize if within buffer distance

    % 4. Time cost
    v_n                 = vecnorm(v, 2,2);
    T                   = sum(1./v_n) * dt;

    % 5. Jerk Cost
    j_mag_sq            = vecnorm(j, 2, 2).^2;         % (N-3) x 1
    v_mag_j             = v_n(1:end-2);                % align size with j
    J                   = sum(j_mag_sq ./ v_mag_j) * dt;      % scalar            

    total_cost = Par.w(1)*L + Par.w(2)*K + Par.w(3)*D + Par.w(4)*T + Par.w(5)*J;

    % Return details for logging or plotting
    global details

    details.length      = L;
    details.curvature   = K;
    details.safety      = D;
    details.time        = T;
    details.jerk        = J;
end