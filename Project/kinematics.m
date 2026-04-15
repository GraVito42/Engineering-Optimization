function [v, a, j, k, dt, dx] = kinematics(P, Par)
%KINEMATICS  Compute kinematic quantities along a path with physical scaling.
%
%   P   : [N x 2] path points in PARAMETRIC coordinates (A→B = unit vector)
%   Par : parameter struct, must contain:
%           Par.dc              - parametric step used to build P
%           Par.LengthReference - real physical distance A→B [m]
%           Par.v_avg           - average cruise speed [m/s]
%
%   Outputs are all in SI units (m, m/s, m/s², m/s³, 1/m, s).

    % --- 1. Parametric arc length (in the normalised coordinate frame) ---
    dP          = diff(P);                          % [N-1 x 2]
    ds_param    = vecnorm(dP, 2, 2);                % parametric step lengths
    L_param     = sum(ds_param);                    % total parametric arc length

    % --- 2. Physical scale factor ---
    % The parametric curve is always 1 (unit vector) [param-unit].
    % Par.LengthReference gives its true physical length reference [m].
    % Every parametric length unit therefore corresponds to lambda metres.
    lambda      = Par.LengthReference/(norm(Par.B - Par.A , 2));         % [m / param-unit]
    dx          = lambda * ds_param;           % physical step between samples [m]

    % --- 3. Physical arc length and travel time ---
    L_phys      = L_param * lambda;                % true path length [m]
    T_phys      = L_phys  / Par.v_avg;             % travel time [s]

    % --- 4. Physical time step between consecutive samples ---
    % For the finite-difference derivatives we need a uniform dt.
    dt          = T_phys * Par.dc;                 % [s]  (uniform step)

    % --- 5. Finite-difference derivatives (physical units) ---
    % Positions are still in parametric coords → scale by lambda first.
    P_phys      = P * lambda;                      % [N x 2]  [m]

    v           = diff(P_phys) / dt;               % [N-1 x 2]  [m/s]
    a           = diff(v)      / dt;               % [N-2 x 2]  [m/s²]
    j           = diff(a)      / dt;               % [N-3 x 2]  [m/s³]

    % --- 6. Curvature [1/m] ---
    v_mid       = (v(1:end-1,:) + v(2:end,:)) / 2;       % midpoint velocity
    cross_prod  = v_mid(:,1).*a(:,2) - v_mid(:,2).*a(:,1);
    k           = abs(cross_prod) ./ (vecnorm(v_mid,2,2).^3 + 1e-6);  % [1/m]

end
