function [g, d] = obstacle_distance(P, Par)
    buffer = Par.d_safe;
    
    % Check if obstacles are defined as a cell array (generic shapes) or numeric array (circles)
    if iscell(Par.obs)
        % --- NEW VERSION: Generic-Shaped Obstacles ---
        num_obs = length(Par.obs);
        g = zeros(num_obs, 1);
        d = zeros(num_obs, 1);
        
        for i = 1:num_obs
            obs_curve = Par.obs{i}; % [Kx2] coordinates of the curve
            
            % Preallocate array to store min distance from each path point to the curve
            min_dist_to_curve = zeros(size(P, 1), 1);
            
            % Compute shortest distance from each point in P to the current obstacle curve
            for j = 1:size(P, 1)
                dist_sq = (obs_curve(:,1) - P(j,1)).^2 + (obs_curve(:,2) - P(j,2)).^2;
                min_dist_to_curve(j) = sqrt(min(dist_sq));
            end
            
            % Check if any path points fall INSIDE the generic polygon
            % (Crucial so that points inside the obstacle aren't seen as "safe" due to positive distance to boundary)
            in = inpolygon(P(:,1), P(:,2), obs_curve(:,1), obs_curve(:,2));
            
            % Treat points inside the obstacle as negative distance to heavily penalize the constraint
            min_dist_to_curve(in) = -min_dist_to_curve(in);
            
            % Overall minimum distance from the path to this obstacle
            min_dist = min(min_dist_to_curve);
            
            % Constraint: buffer - min_dist <= 0 (Negative means safe)
            g(i) = (buffer - min_dist) / buffer;
            
            % For generic shapes, d(i) represents the closest distance to the boundary
            d(i) = min_dist; 
        end
        
    elseif isnumeric(Par.obs)
        % --- Rounded-Shape Obstacles ---
        % centers: [Mx2], radii: [Mx1], P: [Nx2]
        centers = Par.obs(:,1:2);
        radii   = Par.obs(:,3);
        num_obs = size(centers, 1);
        
        g = zeros(num_obs, 1);
        d = zeros(num_obs, 1);
        
        for i = 1:num_obs
            % Distance from every path point to obstacle center
            dist = sqrt((P(:,1) - centers(i,1)).^2 + (P(:,2) - centers(i,2)).^2);
        
            % Surface clearance: subtract radius to get distance to surface
            % Negative if path penetrates the obstacle
            clearance = min(dist) - radii(i);          % [m] surface distance
        
            % Normalized constraint: g=0 at safety boundary, g=1 at surface
            g(i) = (buffer - clearance) / buffer;
        
            % Surface clearance for cost function safety barrier
            d(i) = clearance;                          % [m] NOT center distance
        end
        
    else
        error('Par.obs must be a numeric array for circular obstacles or a cell array for generic shapes.');
    end
end
