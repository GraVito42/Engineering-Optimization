function g = obstacle_distance(centers, radii, P)
    % centers: [Mx2], radii: [Mx1], P: [Nx2]
    num_obs = size(centers, 1);
    g = zeros(num_obs, 1);
    
    for i = 1:num_obs
        % Distance from every point on path to this specific obstacle center
        dist_sq = (P(:,1) - centers(i,1)).^2 + (P(:,2) - centers(i,2)).^2;
        dist = sqrt(dist_sq);
        
        % Constraint: radius - min_dist <= 0  (Negative means safe)
        % fmincon expects g <= 0
        g(i) = radii(i) - min(dist);
    end
end