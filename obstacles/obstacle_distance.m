function d_min = obstacle_distance(center, radius, P)
% OBSTACLE_DISTANCE Computes the clearance from a circular obstacle
% d_min > 0: safe, d_min < 0: collision
    % Euclidean distance from center [xo, yo] to each point in P [Nx2]
    dist_to_center = sqrt((P(:,1) - center(1)).^2 + (P(:,2) - center(2)).^2);
    d_min = min(dist_to_center) - radius;
end