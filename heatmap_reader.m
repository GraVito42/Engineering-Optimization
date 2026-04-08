% =========================================================
% Heatmap Color Contour Extractor - MATLAB Version
% =========================================================

%% --- OPTIONAL: Generate a synthetic color heatmap ---
% Uncomment this section to create 'heatmap_color.png' for testing.
% create_color_heatmap();

%% --- MAIN EXECUTION ---
image_path = 'heatmap_multiple.png';

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