%% =========================================================
%% Island Land vs Sea Heatmap + Contour Extractor
%% =========================================================

clear; clc; close all;

image_path = 'Aruba.png';  % <-- change to your image

%% STEP 1+2: Land/Sea heatmap analysis
[heatmap, blended, land_mask, sea_mask, img_orig] = land_sea_heatmap(image_path);

%% STEP 3: Extract contours from the heatmap's hot zones (land)
[contours, ~, heat_mask] = extract_heat_contours(heatmap);

%% STEP 4: Display everything in ONE window
display_analysis(img_orig, land_mask, sea_mask, heatmap, blended, contours, heat_mask);


%% =========================================================
%%  FUNCTION: land_sea_heatmap
%% =========================================================
function [heatmap, blended, land_mask, sea_mask, img] = land_sea_heatmap(image_path)
    img        = imread(image_path);
    img_double = double(img) / 255;

    hsv = rgb2hsv(img_double);
    H   = hsv(:,:,1);
    S   = hsv(:,:,2);
    V   = hsv(:,:,3);

    sea_mask   = (H >= 0.50 & H <= 0.75) & (S >= 0.15) & (V >= 0.15);
    green_land = (H >= 0.20 & H <= 0.45) & (S >= 0.10);
    sandy_land = (S < 0.25)              & (V > 0.35) & ~sea_mask;
    brown_land = (H >= 0.05 & H <= 0.12) & (S >= 0.15);
    land_mask  = green_land | sandy_land | brown_land;

    se        = strel('disk', 5);
    sea_mask  = imclose(imopen(sea_mask,  se), se);
    land_mask = imclose(imopen(land_mask, se), se);
    land_mask = land_mask & ~sea_mask;
    unknown_mask = ~sea_mask & ~land_mask;

    intensity               = zeros(size(H));
    intensity(land_mask)    = 1.0;
    intensity(sea_mask)     = 0.0;
    intensity(unknown_mask) = 0.45;

    intensity_smooth = imgaussfilt(intensity, 8);
    intensity_smooth = intensity_smooth - min(intensity_smooth(:));
    intensity_smooth = intensity_smooth / max(intensity_smooth(:));

    cmap    = jet(256);
    idx     = gray2ind(intensity_smooth, 256);
    heatmap = ind2rgb(idx, cmap);

    alpha   = 0.65;
    blended = max(0, min(1, alpha * heatmap + (1-alpha) * img_double));

    [~, name] = fileparts(image_path);
    imwrite(uint8(heatmap*255), sprintf('%s_heatmap.png', name));
    imwrite(uint8(blended*255), sprintf('%s_blended.png', name));
end


%% =========================================================
%%  FUNCTION: extract_heat_contours
%% =========================================================
function [contours_out, heatmap_rgb, heat_mask] = extract_heat_contours(heatmap)
    hsv = rgb2hsv(heatmap);
    H   = hsv(:,:,1);
    S   = hsv(:,:,2);
    V   = hsv(:,:,3);

    heat_mask   = ((H <= 0.15) | (H >= 0.85)) & (S >= 0.4) & (V >= 0.4);
    yellow_zone = (H >= 0.12 & H <= 0.20)     & (S >= 0.5) & (V >= 0.5);
    heat_mask   = heat_mask | yellow_zone;

    se        = strel('disk', 4);
    heat_mask = imclose(imopen(heat_mask, se), se);

    boundaries   = bwboundaries(heat_mask, 'noholes');
    contours_out = cell(size(boundaries));
    for k = 1:length(boundaries)
        b = boundaries{k};
        contours_out{k} = [b(:,2), b(:,1)];
    end

    heatmap_rgb = heatmap;
    fprintf('[extract_heat_contours] Found %d heat zone contours.\n', length(contours_out));
end


%% =========================================================
%%  FUNCTION: display_analysis  — ALL PANELS IN ONE WINDOW
%% =========================================================
function display_analysis(img_orig, land_mask, sea_mask, heatmap, blended, contours, ~)

    figure('Name', 'Island Heatmap Analysis', 'Position', [100, 100, 1100, 600]);

    %% Row 1
    subplot(2, 3, 1);
    imshow(img_orig);
    title('Original Image', 'FontSize', 10);
    axis off;

    subplot(2, 3, 2);
    imshow(land_mask);
    title('Land Mask', 'FontSize', 10);
    axis off;

    subplot(2, 3, 3);
    imshow(sea_mask);
    title('Sea Mask', 'FontSize', 10);
    axis off;

    %% Row 2
    subplot(2, 3, 4);
    imshow(heatmap);
    title('Heatmap (Red=Land, Blue=Sea)', 'FontSize', 10);
    axis off;

    subplot(2, 3, 5);
    imshow(blended);
    title('Blended Overlay', 'FontSize', 10);
    axis off;

    subplot(2, 3, 6);
    hold on; grid on; box on;
    axis equal;
    img_h = size(img_orig, 1);
    img_w = size(img_orig, 2);
    xlim([0, img_w]);
    ylim([0, img_h]);
    set(gca, 'YDir', 'reverse');

    colors      = lines(max(length(contours), 1));
    valid_count = 0;
    for i = 1:length(contours)
        c = contours{i};
        if size(c, 1) < 20, continue; end
        valid_count = valid_count + 1;
        plot([c(:,1); c(1,1)], [c(:,2); c(1,2)], ...
             'Color', colors(mod(i-1, size(colors,1))+1, :), ...
             'LineWidth', 1.5);
    end
    xlabel('X (px)', 'FontSize', 9);
    ylabel('Y (px)', 'FontSize', 9);
    title(sprintf('Cartesian Contours (%d zones)', valid_count), 'FontSize', 10);
    hold off;

    sgtitle('Island Land vs Sea Heatmap Analysis', 'FontSize', 13, 'FontWeight', 'bold');
end