% =========================================================
% Heatmap Color Contour Extractor - MATLAB Version
% =========================================================

%% --- OPTIONAL: Generate a synthetic color heatmap ---
% Uncomment this section to create 'heatmap_color.png' for testing.
% create_color_heatmap();

%% --- MAIN EXECUTION ---
image_path = 'heatmap_multiple.png';

% Single image
%land_sea_heatmap(image_path);

% Batch process a folder
imgs = dir('maps/*.jpg');
for i = 1:length(imgs)
    land_sea_heatmap(fullfile('maps', imgs(i).name));
end

[curves, color_img, mask] = extract_color_contours(image_path);

%% --- PLOTTING ---
figure('Position', [100, 100, 1400, 450]);

% Plot 1: Source COLOR Heatmap
subplot(1, 3, 1);
imshow(color_img);
title('Source COLOR Heatmap');
axis off;

% Plot 2: Obstacle Mask
subplot(1, 3, 2);
imshow(mask);
title('Obstacle Mask (Yellow-Red zones)');
axis off;

% Plot 3: Detected Contours overlaid on source
subplot(1, 3, 3);
imshow(color_img);
hold on;
for i = 1:length(curves)
    contour = curves{i};         % [N x 2] array of (x, y) points
    % Close the contour by appending the first point
    x = [contour(:,1); contour(1,1)];
    y = [contour(:,2); contour(1,2)];
    plot(x, y, 'Color', 'cyan', 'LineWidth', 2.5);
end
title(sprintf('Detected Borders (%d obstacles)', length(curves)));
axis off;
hold off;

%% =========================================================
%%  FUNCTION: extract_color_contours
%% =========================================================
function [contours_out, img_rgb, binary_mask] = extract_color_contours(image_path)
    % 1. Load the image (MATLAB reads as RGB)
    img_rgb = imread(image_path);
    if isempty(img_rgb)
        error('Image not found at: %s', image_path);
    end

    % 2. Convert RGB -> HSV
    %    MATLAB's rgb2hsv returns H in [0,1], S in [0,1], V in [0,1]
    hsv = rgb2hsv(img_rgb);
    H = hsv(:,:,1);
    S = hsv(:,:,2);
    V = hsv(:,:,3);

    % 3. Define the heat color range (Yellow-Orange-Red).
    %    OpenCV HSV:  H in [0,180], S in [0,255], V in [0,255]
    %    MATLAB HSV:  H in [0,1],   S in [0,1],   V in [0,1]
    %
    %    OpenCV lower = [0, 150, 50]   -> MATLAB [0/180,  150/255, 50/255]
    %    OpenCV upper = [20, 255, 255] -> MATLAB [20/180, 255/255, 255/255]

    lower_H = 0   / 180;   upper_H = 20  / 180;
    lower_S = 150 / 255;   upper_S = 1.0;
    lower_V = 50  / 255;   upper_V = 1.0;

    % 4. Build the binary mask
    binary_mask = (H >= lower_H & H <= upper_H) & ...
                  (S >= lower_S & S <= upper_S) & ...
                  (V >= lower_V & V <= upper_V);

    % 5. Extract contours using bwboundaries (equivalent to cv2.findContours EXTERNAL)
    %    bwboundaries returns boundaries as {N x 2} cell array of [row, col] = [y, x]
    boundaries = bwboundaries(binary_mask, 'noholes');

    % 6. Reformat to match OpenCV convention: cell array of [x, y] columns
    contours_out = cell(size(boundaries));
    for k = 1:length(boundaries)
        b = boundaries{k};          % [N x 2]: col1=row(y), col2=col(x)
        contours_out{k} = [b(:,2), b(:,1)];  % -> [x, y]
    end
end

%% =========================================================
%%  FUNCTION: create_color_heatmap  (synthetic test image)
%% =========================================================
function create_color_heatmap()
    size_px = 500;
    img = zeros(size_px, size_px, 'uint8');

    % Draw filled circles
    img = draw_filled_circle(img, 150, 150, 80,  200);
    img = draw_filled_circle(img, 300, 350, 100, 255);  % (row,col) = (y,x)
    img = draw_filled_circle(img, 400, 100, 50,  150);

    % Gaussian blur
    img_double = double(img) / 255;
    sigma = 15;
    img_blur = imgaussfilt(img_double, sigma);
    img_blur = uint8(img_blur * 255);

    % Apply JET colormap
    jet_map = colormap(jet(256));
    color_img = ind2rgb(img_blur, jet_map);
    color_img = uint8(color_img * 255);

    imwrite(color_img, 'heatmap_color.png');
    fprintf('Saved heatmap_color.png\n');
end

%% =========================================================
%%  HELPER: draw_filled_circle
%% =========================================================
function img = draw_filled_circle(img, cx, cy, radius, value)
    % cx, cy are (col, row) = (x, y) in image coordinates
    [rows, cols] = size(img);
    [C, R] = meshgrid(1:cols, 1:rows);
    mask = (C - cx).^2 + (R - cy).^2 <= radius^2;
    img(mask) = value;
end

%% =========================================================
%% Land vs Sea Heatmap Generator
%% =========================================================

function land_sea_heatmap(image_path)
    %% 1. Load image
    img = imread(image_path);
    img_double = double(img) / 255;

    %% 2. Convert to HSV for better color segmentation
    hsv = rgb2hsv(img_double);
    H = hsv(:,:,1);
    S = hsv(:,:,2);
    V = hsv(:,:,3);

    %% 3. Detect SEA (blue tones)
    % Hue around 0.55-0.75 (cyan to blue), decent saturation
    sea_mask = (H >= 0.50 & H <= 0.75) & ...
               (S >= 0.15) & ...
               (V >= 0.15);

    %% 4. Detect LAND (greens, yellows, browns, sandy tones)
    % Green land
    green_land = (H >= 0.20 & H <= 0.45) & (S >= 0.10);
    % Sandy/desert/urban (low saturation, mid-high brightness)
    sandy_land = (S < 0.25) & (V > 0.35) & ~sea_mask;
    % Brown/earthy tones
    brown_land = (H >= 0.05 & H <= 0.12) & (S >= 0.15);

    land_mask = green_land | sandy_land | brown_land;

    %% 5. Clean up masks with morphological ops
    se = strel('disk', 5);
    sea_mask  = imopen(sea_mask,  se);
    sea_mask  = imclose(sea_mask, se);
    land_mask = imopen(land_mask,  se);
    land_mask = imclose(land_mask, se);

    % Resolve overlap: sea takes priority if both detected
    land_mask = land_mask & ~sea_mask;

    % Unknown/unclassified pixels
    unknown_mask = ~sea_mask & ~land_mask;

    %% 6. Build intensity map for heatmap
    % Land = high intensity (hot), Sea = low intensity (cold), Unknown = mid
    intensity = zeros(size(H));
    intensity(land_mask)    = 1.0;   % land   → hot (red)
    intensity(sea_mask)     = 0.0;   % sea    → cold (blue)
    intensity(unknown_mask) = 0.45;  % unknown → neutral

    %% 7. Smooth the intensity for a natural heatmap look
    intensity_smooth = imgaussfilt(intensity, 8);
    intensity_smooth = intensity_smooth - min(intensity_smooth(:));
    intensity_smooth = intensity_smooth / max(intensity_smooth(:));

    %% 8. Apply JET colormap
    cmap    = jet(256);
    idx     = gray2ind(intensity_smooth, 256);
    heatmap = ind2rgb(idx, cmap);

    %% 9. Blend heatmap over original
    alpha   = 0.65;
    blended = alpha * heatmap + (1 - alpha) * img_double;
    blended = max(0, min(1, blended));  % clamp

    %% 10. Save outputs
    [~, name, ~] = fileparts(image_path);
    imwrite(uint8(heatmap * 255), sprintf('%s_heatmap.png',  name));
    imwrite(uint8(blended * 255), sprintf('%s_blended.png',  name));
    imwrite(uint8(land_mask) * 255, sprintf('%s_land_mask.png', name));
    imwrite(uint8(sea_mask)  * 255, sprintf('%s_sea_mask.png',  name));
    fprintf('Saved outputs for: %s\n', image_path);

    %% 11. Display
    figure('Position', [50, 50, 1600, 800]);

    subplot(2, 3, 1);
    imshow(img);
    title('Original Map');
    axis off;

    subplot(2, 3, 2);
    imshow(land_mask);
    title('Land Mask');
    axis off;

    subplot(2, 3, 3);
    imshow(sea_mask);
    title('Sea Mask');
    axis off;

    subplot(2, 3, 4);
    imshow(heatmap);
    title('Heatmap (Red=Land, Blue=Sea)');
    axis off;

    subplot(2, 3, 5);
    imshow(blended);
    title('Blended (Heatmap + Original)');
    axis off;

    subplot(2, 3, 6);
    % Coverage stats bar chart
    total  = numel(H);
    land_pct    = 100 * sum(land_mask(:))    / total;
    sea_pct     = 100 * sum(sea_mask(:))     / total;
    unknown_pct = 100 * sum(unknown_mask(:)) / total;
    bar([land_pct, sea_pct, unknown_pct], 'FaceColor', 'flat', ...
        'CData', [0.8 0.3 0.1; 0.1 0.4 0.9; 0.6 0.6 0.6]);
    set(gca, 'XTickLabel', {'Land', 'Sea', 'Unknown'});
    ylabel('Coverage (%)');
    title('Detection Stats');
    ylim([0 100]);
    grid on;

    sgtitle(sprintf('Land vs Sea Analysis: %s', image_path), 'FontSize', 14);
end