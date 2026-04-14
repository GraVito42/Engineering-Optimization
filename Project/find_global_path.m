function [raw_path, map] = find_global_path(Par)
% FIND_GLOBAL_PATH Creates an occupancy map and finds a safe corridor using A*
% after it created an occupancy map matrix, adujusting the obs with the buffer
% A* is an algorithm that find the sortest path to connet two matrix cells
% avoiding the obstacles

    fprintf('Building occupancy map and running A*...\n');
    
    % Find boundaries
    nodes_x = [Par.A(1), Par.B(1)];
    nodes_y = [Par.A(2), Par.B(2)];
    for i = 1:length(Par.obs)
        if ~isempty(Par.obs{i})
            nodes_x = [nodes_x, Par.obs{i}(:,1)'];
            nodes_y = [nodes_y, Par.obs{i}(:,2)'];
        end
    end
    
    buffer = 50;                % needs to be adjusted to the physical lenght over resolution 
    min_x = floor(min(nodes_x)) - buffer;  
    max_x = ceil(max(nodes_x)) + buffer;

    min_y = floor(min(nodes_y)) - buffer;  
    max_y = ceil(max(nodes_y)) + buffer;

    mapWidth  = max_x - min_x;
    mapHeight = max_y - min_y;
    
    % Create map
    map = binaryOccupancyMap(mapWidth, mapHeight, 1);
    map.LocalOriginInWorld = [min_x, min_y]; 

    % Insert and inflate obstacles
    for i = 1:length(Par.obs)
        if ~isempty(Par.obs{i})
            setOccupancy(map, Par.obs{i}, 1);
        end
    end
    inflate(map, Par.d_safe);    % applies our buffer to the matrix

    % Run A*
    planner = plannerAStarGrid(map);

    start_grid = world2grid(map, Par.A(:)');
    goal_grid  = world2grid(map, Par.B(:)');

    path_grid = plan(planner, start_grid, goal_grid);

    if isempty(path_grid)
        raw_path = NaN;
    else
        raw_path = grid2world(map, path_grid);
    end

    
end