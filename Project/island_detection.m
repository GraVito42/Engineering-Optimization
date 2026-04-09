function [island_curves, img_orig, heatmap, blended] = island_detection(image_path, show_plots)
%ISLAND_DETECTION  Detect land/island contours from a map image.
%
%   [island_curves, img_orig, heatmap, blended] = island_detection(image_path, show_plots)
%
%   INPUTS:
%     image_path  - (string) path to the input map image (e.g. 'Aruba.png')
%     show_plots  - (logical) true = display the 6-panel analysis figure, false = silent
%
%   OUTPUTS:
%     island_curves - (cell array) each cell contains an [Kx2] array of [x, y]
%                     pixel coordinates forming one island/land contour.
%                     Compatible with the generic-shape branch of obstacle_distance().
%     img_orig      - (HxWx3 uint8)  original image
%     heatmap       - (HxWx3 double) jet-coloured land/sea heatmap  [0,1]
%     blended       - (HxWx3 double) heatmap blended over original  [0,1]

    if nargin < 2
        show_plots = true;
    end

    %% --- Stage 1 & 2: build land/sea masks and heatmap ---
    [heatmap, blended, land_mask, sea_mask, img_orig] = build_heatmap(image_path);

    %% --- Stage 3: extract island contours from the heatmap hot zones ---
    [island_curves, heat_mask] = extract_island_contours(heatmap);

    %% --- Stage 4 (optional): visualise ---
    if show_plots
        display_analysis(img_orig, land_mask, sea_mask, ...
                         heatmap, blended, island_curves, heat_mask);
    end
end


%% =========================================================
%%  LOCAL FUNCTION: build_heatmap
%% =========================================================
function [heatmap, blended, land_mask, sea_mask, img] = build_heatmap(image_path)
    img        = imread(image_path);
    img_double = double(img) / 255;

    hsv = rgb2hsv(img_double);
    H   = hsv(:,:,1);
    S   = hsv(:,:,2);
    V   = hsv(:,:,3);

    % --- colour-based segmentation ---
    sea_mask   = (H >= 0.50 & H <= 0.75) & (S >= 0.15) & (V >= 0.15);
    green_land = (H >= 0.20 & H <= 0.45) & (S >= 0.10);
    sandy_land = (S < 0.25)              & (V > 0.35) & ~sea_mask;
    brown_land = (H >= 0.05 & H <= 0.12) & (S >= 0.15);
    land_mask  = green_land | sandy_land | brown_land;

    % --- morphological clean-up ---
    se        = strel('disk', 5);
    sea_mask  = imclose(imopen(sea_mask,  se), se);
    land_mask = imclose(imopen(land_mask, se), se);
    land_mask = land_mask & ~sea_mask;
    unknown_mask = ~sea_mask & ~land_mask;

    % --- intensity map: 1=land, 0=sea, 0.45=unknown ---
    intensity               = zeros(size(H));
    intensity(land_mask)    = 1.0;
    intensity(sea_mask)     = 0.0;
    intensity(unknown_mask) = 0.45;

    intensity_smooth = imgaussfilt(intensity, 8);
    intensity_smooth = intensity_smooth - min(intensity_smooth(:));
    intensity_smooth = intensity_smooth / max(intensity_smooth(:));

    % --- map to jet colourmap ---
    cmap    = jet(256);
    idx     = gray2ind(intensity_smooth, 256);
    heatmap = ind2rgb(idx, cmap);

    alpha   = 0.65;
    blended = max(0, min(1, alpha * heatmap + (1-alpha) * img_double));

    % --- save side-products ---
    [~, name] = fileparts(image_path);
    imwrite(uint8(heatmap * 255), sprintf('%s_heatmap.png', name));
    imwrite(uint8(blended  * 255), sprintf('%s_blended.png',  name));
end


%% =========================================================
%%  LOCAL FUNCTION: extract_island_contours
%% =========================================================
function [island_curves, heat_mask] = extract_island_contours(heatmap)
%   Returns a cell array where each cell is one closed land contour [Kx2]
%   in the format expected by obstacle_distance() (generic-shape branch):
%       obs_curve(:,1) = x  (column pixel)
%       obs_curve(:,2) = y  (row    pixel)

    hsv = rgb2hsv(heatmap);
    H   = hsv(:,:,1);
    S   = hsv(:,:,2);
    V   = hsv(:,:,3);

    % hot colours in the jet map correspond to land
    heat_mask   = ((H <= 0.15) | (H >= 0.85)) & (S >= 0.4) & (V >= 0.4);
    yellow_zone = (H >= 0.12 & H <= 0.20)     & (S >= 0.5) & (V >= 0.5);
    heat_mask   = heat_mask | yellow_zone;

    se        = strel('disk', 4);
    heat_mask = imclose(imopen(heat_mask, se), se);

    boundaries   = bwboundaries(heat_mask, 'noholes');
    island_curves = cell(size(boundaries));
    for k = 1:length(boundaries)
        b = boundaries{k};
        % Convert from [row, col] (bwboundaries convention) to [x, y] = [col, row]
        island_curves{k} = [b(:,2), -b(:,1)];
    end

    fprintf('[island_detection] Found %d island/land contours.\n', length(island_curves));
end


%% =========================================================
%%  LOCAL FUNCTION: display_analysis
%% =========================================================
function display_analysis(img_orig, land_mask, sea_mask, heatmap, blended, island_curves, heat_mask)

    figure('Name', 'Island Heatmap Analysis', 'Position', [100, 100, 1100, 600]);

    subplot(2, 3, 1);
    imshow(img_orig);
    title('Original Image', 'FontSize', 10);

    subplot(2, 3, 2);
    imshow(land_mask);
    title('Land Mask', 'FontSize', 10);

    subplot(2, 3, 3);
    imshow(sea_mask);
    title('Sea Mask', 'FontSize', 10);

    subplot(2, 3, 4);
    imshow(heatmap);
    title('Heatmap  (Red = Land, Blue = Sea)', 'FontSize', 10);

    subplot(2, 3, 5);
    imshow(blended);
    title('Blended Overlay', 'FontSize', 10);

    subplot(2, 3, 6);
    hold on; grid on; box on; axis equal;
    img_h = size(img_orig, 1);
    img_w = size(img_orig, 2);
    xlim([0, img_w]);  ylim([0, img_h]);
    % set(gca, 'YDir', 'reverse');

    colors      = lines(max(length(island_curves), 1));
    valid_count = 0;
    for i = 1:length(island_curves)
        c = island_curves{i};
        if size(c, 1) < 20, continue; end
        valid_count = valid_count + 1;
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], ...
             'Color', colors(mod(i-1, size(colors,1))+1, :), ...
             'LineWidth', 1.5);
    end
    xlabel('X (px)');  ylabel('Y (px)');
    title(sprintf('Island Contours  (%d zones)', valid_count), 'FontSize', 10);
    hold off;

    sgtitle('Island Land vs Sea Heatmap Analysis', 'FontSize', 13, 'FontWeight', 'bold');
end