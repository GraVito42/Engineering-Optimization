function [total_cost, details] = bernstein_cost(P, dt)
% BERNSTEIN_COST Computes the analytical cost of a 2D trajectory
%
% INPUTS:
%   P  : [Nx2] Matrix of trajectory points (x, y)
%   dt : Scalar, time step between samples
%
% OUTPUTS:
%   total_cost : Weighted sum of length and energy
%   details    : Structure containing individual cost components

    % 1. Compute Derivatives (Velocity and Acceleration)
    % We use central differences for better accuracy or simple diff
    v = diff(P) / dt;       % Velocity vectors [ (N-1) x 2 ]
    a = diff(v) / dt;       % Acceleration vectors [ (N-2) x 2 ]

    % 2. Length Cost (L2 Norm of velocity)
    % Sum of infinitesimal displacements: ds = |v| * dt
    segment_lengths = sqrt(sum(v.^2, 2));
    L_cost = sum(segment_lengths) * dt;

    % 3. Smoothness/Energy Cost (Squared L2 Norm of acceleration)
    % E = integral( |a|^2 * dt )
    accel_squared = sum(a.^2, 2);
    E_cost = sum(accel_squared) * dt;

    % 4. Weighted Total Cost
    % w1 and w2 can be tuned: w1 for urgency, w2 for battery/smoothness
    w1 = 1.0; 
    w2 = 0.5; 
    total_cost = w1 * L_cost + w2 * E_cost;

    % Return details for logging or plotting
    details.length = L_cost;
    details.energy = E_cost;
end