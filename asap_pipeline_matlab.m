%% ASAP burn severity analysis pipeline
% Integrates the three ASAP preprocessing steps into a single script
%
%   Step 1. Definition of the analysis area
%   Step 2. Preparation of analysis imagery
%   Step 3. Phenology detrending
%
% Outputs from each step are saved separately to:
%   01_ASAP_step1, 02_ASAP_step2, and 03_ASAP_step3
%
% Notes
% - By default, the project root is set to the folder containing this
%   script. If needed, manually edit "path_root" below.

clear; close all; clc;

%% ------------------------------------------------------------------------
% 1. User settings
% -------------------------------------------------------------------------

% Project root
% If this script is placed in the project root, the default below is enough.
% Otherwise, set "path_root" manually.
path_root = fileparts(mfilename('fullpath'));
% path_root = 'T:\tjsung\Research\L8_BurnSeverity_CONUS\Github\';
% path_root = '/share/wildfire-2/tjsung/Research/L8_BurnSeverity_CONUS/Github/';

% Main folders
paths = get_project_paths(path_root);

% Processing options
cfg.clear_ratio_threshold = 0.70;   % Minimum clear ratio for general Landsat scenes
cfg.af_buffer_size        = 201;    % Buffered AF window size (square kernel)
cfg.roi_padding_pixels    = 100;    % Padding used for fallback ROI
cfg.vr_step_pixels        = 100;    % Expansion size per VR iteration
cfg.vr_max_steps          = 5;      % Maximum VR expansion iterations
cfg.min_lr_samples        = 10;     % Minimum unique samples for linear fitting
cfg.reflectance_min_dn    = 7273;   % Valid Landsat DN lower bound
cfg.reflectance_max_dn    = 43636;  % Valid Landsat DN upper bound
cfg.firms_window_days     = 28;     % Time window after ignition for FIRMS detections

% Landsat constants
cfg.qa_clear_value        = 21824;
cfg.qa_exclude_value      = 21952;
cfg.reflectance_scale     = 0.0000275;
cfg.reflectance_offset    = -0.2;

% Figure settings
cfg.figure_position       = [-1000, 0, 1000, 1000];
cfg.map_font_size         = 25;
cfg.colorbar_font_size    = 18;
cfg.title_font_size       = 28;
cfg.step1_title_font_size = 24;
cfg.fire_map_bg_color     = [0.18, 0.28, 0.18];
cfg.si_clim               = [0.0, 0.7];
cfg.dsi_clim              = [-0.7, 0.7];

% Index names
cfg.list_si   = {'NDVI', 'NBR'};
cfg.list_dsi  = {'dNDVI', 'dNBR', 'rdNBR', 'RBR'};
cfg.list_pdsi = {'pdNDVI', 'pdNBR', 'prdNBR', 'pRBR'};

% Create output folders
create_output_dirs(paths);

%% ------------------------------------------------------------------------
% 2. Read common inputs
% -------------------------------------------------------------------------

metadata = readtable(fullfile(paths.data, 'metadata.xlsx'), 'PreserveVariableNames', 1);

% FIRMS active fire products used in the original step 1 workflow
%   - MODIS: Terra/Aqua
%   - VIIRS: S-NPP/JPSS-1/JPSS-2
firms_modis = readtable(fullfile(paths.data, 'FIRMS', 'fire_M-C61.csv'), ...
    'PreserveVariableNames', 1);
firms_viirs = readtable(fullfile(paths.data, 'FIRMS', 'fire_V-C2.csv'), ...
    'PreserveVariableNames', 1);

%% ------------------------------------------------------------------------
% 3. Main loop: process each wildfire case
% -------------------------------------------------------------------------

n_case = height(metadata);

for ii = 1:n_case

    case_name = metadata.Incid_Name{ii};
    ig_date   = metadata.Ig_Date(ii);
    post_date = metadata.Post_Date(ii);
    pre_date  = metadata{ii, 'Pre_Date(const)'};

    fprintf('\n[%d/%d] Processing case: %s\n', ii, n_case, case_name);

    try
        %% =================================================================
        % STEP 1. Definition of the analysis area
        % ==================================================================

        [MTBS_BNDY, LON, LAT] = read_mtbs_boundary(case_name, paths.data);
        NLCD_LC = read_nlcd_reclassified(case_name, paths.data);

        [MODIS_AF, pxid_modis] = rasterize_firms_detections( ...
            firms_modis, ig_date, LON, LAT, NLCD_LC, cfg);
        [VIIRS_AF, pxid_viirs] = rasterize_firms_detections( ...
            firms_viirs, ig_date, LON, LAT, NLCD_LC, cfg);

        AF_BUFFER = imdilate(MODIS_AF | VIIRS_AF, strel('square', cfg.af_buffer_size));
        af_boundaries = bwboundaries(AF_BUFFER, 'noholes');
        roi_bboxes = boundaries_to_bboxes(af_boundaries);

        target_idx = select_center_roi_bbox(roi_bboxes, size(MTBS_BNDY));

        [ROI, ROI_AF, target_bbox] = define_roi( ...
            MTBS_BNDY, roi_bboxes, target_idx, cfg.roi_padding_pixels);

        VR = define_vr( ...
            ROI, ROI_AF, NLCD_LC, target_bbox, cfg.vr_step_pixels, cfg.vr_max_steps);

        save_fire_map( ...
            fullfile(paths.step1_figures, [case_name, '.png']), ...
            case_name, LON, LAT, MTBS_BNDY, pxid_modis, pxid_viirs, ...
            af_boundaries, mask_to_bbox(ROI), mask_to_bbox(VR), cfg);

        save(fullfile(paths.step1_dataset, [case_name, '.mat']), ...
            'MTBS_BNDY', 'NLCD_LC', 'ROI', 'VR', 'LON', 'LAT', '-v7.3');

        %% =================================================================
        % STEP 2. Preparation of analysis imagery
        % ==================================================================

        crop_bbox = mask_to_bbox(VR);

        [LAT, LON, MTBS_BNDY, NLCD_LC, ROI, VR] = crop_base_layers( ...
            LAT, LON, MTBS_BNDY, NLCD_LC, ROI, VR, crop_bbox);

        burn_boundaries = bwboundaries(MTBS_BNDY);

        [L8_SI, L8_TIME] = prepare_landsat_spectral_indices( ...
            case_name, pre_date, post_date, paths, crop_bbox, ...
            LON, LAT, NLCD_LC, burn_boundaries, cfg);

        save(fullfile(paths.step2_dataset, [case_name, '.mat']), ...
            'L8_SI', 'L8_TIME', 'MTBS_BNDY', 'NLCD_LC', 'ROI', 'VR', 'LON', 'LAT', '-v7.3');

        %% =================================================================
        % STEP 3. Phenology detrending
        % ==================================================================

        [L8_dSI, L8_pdSI, L8_TIME_pre, COEF] = run_phenology_detrending( ...
            case_name, ig_date, post_date, L8_SI, L8_TIME, MTBS_BNDY, ...
            NLCD_LC, VR, LON, LAT, burn_boundaries, paths, cfg);

        save(fullfile(paths.step3_dataset, [case_name, '.mat']), ...
            'L8_dSI', 'L8_pdSI', 'L8_TIME_pre', 'COEF', ...
            'MTBS_BNDY', 'NLCD_LC', 'ROI', 'VR', 'LON', 'LAT', '-v7.3');

        fprintf('  -> Finished successfully: %s\n', case_name);

    catch ME
        warning('Case "%s" failed.\n%s', case_name, getReport(ME, 'basic', 'hyperlinks', 'off'));
    end
end

%% ------------------------------------------------------------------------
% Local functions
% -------------------------------------------------------------------------

function paths = get_project_paths(path_root)
% Build all project paths from the project root.

    paths.root           = path_root;
    paths.data           = fullfile(path_root, '00_Data');

    paths.step1          = fullfile(path_root, '01_ASAP_step1');
    paths.step2          = fullfile(path_root, '02_ASAP_step2');
    paths.step3          = fullfile(path_root, '03_ASAP_step3');

    paths.step1_dataset  = fullfile(paths.step1, 'Dataset');
    paths.step1_figures  = fullfile(paths.step1, 'Figures');
    paths.step2_dataset  = fullfile(paths.step2, 'Dataset');
    paths.step2_figures  = fullfile(paths.step2, 'Figures');
    paths.step3_dataset  = fullfile(paths.step3, 'Dataset');
    paths.step3_figures  = fullfile(paths.step3, 'Figures');
end

function create_output_dirs(paths)
% Create output folders if they do not already exist.

    folder_list = { ...
        paths.step1_dataset, paths.step1_figures, ...
        paths.step2_dataset, paths.step2_figures, ...
        paths.step3_dataset, paths.step3_figures};

    for i = 1:numel(folder_list)
        ensure_dir(folder_list{i});
    end
end

function ensure_dir(folder_path)
% Create a folder only if it does not already exist.

    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end
end

function [MTBS_BNDY, LON, LAT] = read_mtbs_boundary(case_name, path_data)
% Read the MTBS/WFIGS boundary raster and build matching longitude/latitude
% grids from the GeoTIFF spatial reference.

    tif_path = fullfile(path_data, 'MTBS_WFIGS', case_name);

    info = geotiffinfo(tif_path);
    [lon_grid, lat_grid] = pixcenters(info);

    LON = repmat(lon_grid, numel(lat_grid), 1);
    LAT = repmat(lat_grid(:), 1, numel(lon_grid));

    MTBS_BNDY = geotiffread(tif_path);
end

function NLCD_LC = read_nlcd_reclassified(case_name, path_data)
% Read the NLCD raster for one wildfire case and reclassify it into the
% six vegetation classes used in ASAP.
%
% Reclassified classes
%   0: Excluded / non-vegetated
%   1: Deciduous forest
%   2: Evergreen forest
%   3: Mixed forest
%   4: Shrubland
%   5: Herbaceous
%   6: Planted/cultivated

    folder_nlcd = fullfile(path_data, 'NLCD');
    file_info = dir(fullfile(folder_nlcd, [case_name, '_*.tif']));

    if isempty(file_info)
        error('NLCD file not found for case: %s', case_name);
    end

    tif_path = fullfile(folder_nlcd, file_info(1).name);
    NLCD_LC = geotiffread(tif_path);

    NLCD_LC(NLCD_LC < 41 | NLCD_LC > 82) = 0;
    NLCD_LC(NLCD_LC == 41) = 1;
    NLCD_LC(NLCD_LC == 42) = 2;
    NLCD_LC(NLCD_LC == 43) = 3;
    NLCD_LC(NLCD_LC >= 51 & NLCD_LC <= 52) = 4;
    NLCD_LC(NLCD_LC >= 71 & NLCD_LC <= 74) = 5;
    NLCD_LC(NLCD_LC >= 81 & NLCD_LC <= 82) = 6;
end

function [AF_map, pxid] = rasterize_firms_detections(firms_table, ig_date, LON, LAT, NLCD_LC, cfg)
% Convert FIRMS point detections into a raster aligned with the MTBS/NLCD
% grid.
%
% Processing steps
%   1. Keep detections inside the temporal range
%   2. Keep detections inside the target image extent
%   3. Assign each detection to the nearest raster cell
%   4. Remove detections over excluded land-cover classes

    time_mask = firms_table.acq_date >= ig_date & ...
                firms_table.acq_date <= ig_date + days(cfg.firms_window_days);
    firms_sub = firms_table(time_mask, :);

    lat_min = min(LAT(:));
    lat_max = max(LAT(:));
    lon_min = min(LON(:));
    lon_max = max(LON(:));

    firms_sub = firms_sub( ...
        firms_sub.latitude  >= lat_min & firms_sub.latitude  <= lat_max & ...
        firms_sub.longitude >= lon_min & firms_sub.longitude <= lon_max, :);

    AF_map = zeros(size(LON), 'uint8');

    if isempty(firms_sub)
        pxid = [];
        return;
    end

    lon_axis = LON(1, :);
    lat_axis = LAT(:, 1);

    for i = 1:height(firms_sub)
        [~, col_idx] = min(abs(lon_axis - firms_sub.longitude(i)));
        [~, row_idx] = min(abs(lat_axis - firms_sub.latitude(i)));
        AF_map(row_idx, col_idx) = 1;
    end

    AF_map(NLCD_LC == 0) = 0;
    pxid = find(AF_map == 1);
end

function bboxes = boundaries_to_bboxes(boundaries)
% Convert boundary coordinates returned by bwboundaries into rectangular
% bounding boxes [row_min, row_max, col_min, col_max].

    if isempty(boundaries)
        bboxes = zeros(0, 4);
        return;
    end

    bboxes = zeros(numel(boundaries), 4);
    for i = 1:numel(boundaries)
        boundary = boundaries{i};
        bboxes(i, :) = [ ...
            min(boundary(:, 1)), max(boundary(:, 1)), ...
            min(boundary(:, 2)), max(boundary(:, 2))];
    end
end

function target_idx = select_center_roi_bbox(roi_bboxes, img_size)
% Select the ROI candidate whose bounding box contains the image center.
% If no candidate contains the center, return 0 and use the fallback ROI.

    target_idx = 0;

    if isempty(roi_bboxes)
        return;
    end

    center_row = (img_size(1) + 1) / 2;
    center_col = (img_size(2) + 1) / 2;

    for i = 1:size(roi_bboxes, 1)
        row_min = roi_bboxes(i, 1);
        row_max = roi_bboxes(i, 2);
        col_min = roi_bboxes(i, 3);
        col_max = roi_bboxes(i, 4);

        if center_row >= row_min && center_row <= row_max && ...
           center_col >= col_min && center_col <= col_max
            target_idx = i;
            return;
        end
    end
end

function [ROI, ROI_AF, target_bbox] = define_roi(MTBS_BNDY, roi_bboxes, target_idx, roi_padding_pixels)
% Define the rectangular ROI.
%
% If a target AF-based ROI exists, use that bounding box.
% Otherwise, build a fallback ROI from the MTBS burn perimeter with padding.

    ROI    = zeros(size(MTBS_BNDY), 'uint8');
    ROI_AF = zeros(size(MTBS_BNDY), 'uint8');

    if target_idx == 0
        target_bbox = compute_mask_bbox(MTBS_BNDY == 1, roi_padding_pixels);
        if isempty(target_bbox)
            error('Fallback ROI cannot be defined because MTBS_BNDY is empty.');
        end
    else
        target_bbox = roi_bboxes(target_idx, :);
    end

    ROI(target_bbox(1):target_bbox(2), target_bbox(3):target_bbox(4)) = 1;

    for i = 1:size(roi_bboxes, 1)
        bbox = roi_bboxes(i, :);
        ROI_AF(bbox(1):bbox(2), bbox(3):bbox(4)) = 1;
    end
end

function VR = define_vr(ROI, ROI_AF, NLCD_LC, target_bbox, vr_step_pixels, vr_max_steps)
% Define the vicinity region (VR) by progressively expanding a rectangle
% around the target ROI until all vegetation classes in the ROI are also
% represented in the VR.
%
% The VR excludes all AF-based ROI candidates.

    VR = zeros(size(ROI), 'uint8');

    roi_classes = unique(NLCD_LC(ROI == 1));
    roi_classes(roi_classes == 0) = [];

    for j = 1:vr_max_steps
        row_min = max(target_bbox(1) - vr_step_pixels * j, 1);
        row_max = min(target_bbox(2) + vr_step_pixels * j, size(ROI, 1));
        col_min = max(target_bbox(3) - vr_step_pixels * j, 1);
        col_max = min(target_bbox(4) + vr_step_pixels * j, size(ROI, 2));

        VR = zeros(size(ROI), 'uint8');
        VR(row_min:row_max, col_min:col_max) = 1;

        % Exclude AF-based candidate regions and the target ROI itself
        VR(ROI_AF == 1) = 0;
        VR(ROI == 1) = 0;

        vr_classes = unique(NLCD_LC(VR == 1));
        vr_classes(vr_classes == 0) = [];

        if all(ismember(roi_classes, vr_classes))
            break;
        end
    end
end

function bbox = compute_mask_bbox(mask, padding)
% Convert a binary mask into a padded bounding box:
% [row_min, row_max, col_min, col_max].

    [rows, cols] = find(mask);

    if isempty(rows)
        bbox = [];
        return;
    end

    bbox = [ ...
        max(min(rows) - padding, 1), ...
        min(max(rows) + padding, size(mask, 1)), ...
        max(min(cols) - padding, 1), ...
        min(max(cols) + padding, size(mask, 2))];
end

function bbox = mask_to_bbox(mask)
% Convert a binary mask into its tight bounding box:
% [row_min, row_max, col_min, col_max].

    bbox = compute_mask_bbox(mask == 1, 0);

    if isempty(bbox)
        error('Bounding box cannot be generated because the mask is empty.');
    end
end

function [LAT_c, LON_c, MTBS_BNDY_c, NLCD_LC_c, ROI_c, VR_c] = ...
            crop_base_layers(LAT, LON, MTBS_BNDY, NLCD_LC, ROI, VR, crop_bbox)
% Crop all base layers using the same bounding box.

    row_min = crop_bbox(1);
    row_max = crop_bbox(2);
    col_min = crop_bbox(3);
    col_max = crop_bbox(4);

    LAT_c       = LAT(row_min:row_max, col_min:col_max);
    LON_c       = LON(row_min:row_max, col_min:col_max);
    MTBS_BNDY_c = MTBS_BNDY(row_min:row_max, col_min:col_max);
    NLCD_LC_c   = NLCD_LC(row_min:row_max, col_min:col_max);
    ROI_c       = ROI(row_min:row_max, col_min:col_max);
    VR_c        = VR(row_min:row_max, col_min:col_max);
end

function [L8_SI, L8_TIME] = prepare_landsat_spectral_indices( ...
            case_name, pre_date, post_date, paths, crop_bbox, ...
            LON, LAT, NLCD_LC, boundaries, cfg)
% Read Landsat scenes for one case, apply QA and land-cover filtering,
% compute NDVI and NBR, and save a figure for every retained scene.
%
% Retention rule
%   - Always keep the metadata-based pre- and post-fire scenes
%   - For all other dates, keep only scenes whose clear ratio exceeds the
%     specified threshold

    list_l8 = dir(fullfile(paths.data, 'Landsat', [case_name, '*.tif']));
    if isempty(list_l8)
        error('No Landsat files found for case: %s', case_name);
    end

    file_names = string({list_l8.name});
    date_str = extractBetween(file_names, strlength(file_names) - 11, strlength(file_names) - 4);
    date_num = str2double(date_str);
    [~, sort_idx] = sort(date_num);
    list_l8 = list_l8(sort_idx);

    si_list = cell(numel(list_l8), 1);
    time_list = NaT(numel(list_l8), 1);
    keep_count = 0;

    row_min = crop_bbox(1);
    row_max = crop_bbox(2);
    col_min = crop_bbox(3);
    col_max = crop_bbox(4);

    for i = 1:numel(list_l8)
        file_name = list_l8(i).name;
        img_date = datetime(file_name(end-11:end-4), 'InputFormat', 'yyyyMMdd');

        l8_raw = single(geotiffread(fullfile(paths.data, 'Landsat', file_name)));
        l8_raw = l8_raw(row_min:row_max, col_min:col_max, :);

        reflectance_bands = l8_raw(:, :, 1:6);
        qa_band = l8_raw(:, :, end);

        keep_scene = (img_date == post_date) || (img_date == pre_date);
        if ~keep_scene
            clear_ratio = compute_clear_ratio(qa_band, cfg);
            if clear_ratio < cfg.clear_ratio_threshold
                continue;
            end
        end

        reflectance_bands = apply_landsat_qc( ...
            reflectance_bands, qa_band, NLCD_LC, cfg);

        [NDVI, NBR] = calculate_spectral_indices(reflectance_bands);

        keep_count = keep_count + 1;
        si_list{keep_count} = cat(3, NDVI, NBR);
        time_list(keep_count) = img_date;

        save_index_map(LON, LAT, NDVI, boundaries, case_name, img_date, ...
            'NDVI', paths.step2, cfg.si_clim, 'si', cfg);
        save_index_map(LON, LAT, NBR, boundaries, case_name, img_date, ...
            'NBR', paths.step2, cfg.si_clim, 'si', cfg);
    end

    if keep_count == 0
        error('No Landsat scenes were retained for case: %s', case_name);
    end

    L8_SI = cat(4, si_list{1:keep_count});
    L8_TIME = time_list(1:keep_count);
end

function reflectance_bands = apply_landsat_qc(reflectance_bands, qa_band, NLCD_LC, cfg)
% Apply Landsat QA filtering, land-cover masking, and DN-to-reflectance
% conversion.

    qa_mask = qa_band ~= cfg.qa_clear_value;
    lc_mask = NLCD_LC == 0;

    reflectance_bands(repmat(qa_mask, [1, 1, size(reflectance_bands, 3)])) = nan;
    reflectance_bands(repmat(lc_mask, [1, 1, size(reflectance_bands, 3)])) = nan;

    reflectance_bands(reflectance_bands < cfg.reflectance_min_dn | ...
                      reflectance_bands > cfg.reflectance_max_dn) = nan;

    reflectance_bands = reflectance_bands * cfg.reflectance_scale + cfg.reflectance_offset;
end

function [NDVI, NBR] = calculate_spectral_indices(reflectance_bands)
% Calculate NDVI and NBR from Landsat surface reflectance bands.
%
% Landsat band order used here
%   band 3: red
%   band 4: NIR
%   band 6: SWIR2

    NDVI = (reflectance_bands(:, :, 4) - reflectance_bands(:, :, 3)) ./ ...
           (reflectance_bands(:, :, 4) + reflectance_bands(:, :, 3));

    NBR = (reflectance_bands(:, :, 4) - reflectance_bands(:, :, 6)) ./ ...
          (reflectance_bands(:, :, 4) + reflectance_bands(:, :, 6));
end

function clear_ratio = compute_clear_ratio(qa_band, cfg)
% Compute the clear ratio used in the original pipeline.

    numerator = sum(qa_band == cfg.qa_clear_value, 'all');
    denominator = sum(qa_band > 1 & qa_band ~= cfg.qa_exclude_value, 'all');

    if denominator == 0
        clear_ratio = 0;
    else
        clear_ratio = numerator / denominator;
    end
end

function [L8_dSI, L8_pdSI, L8_TIME_pre, COEF] = run_phenology_detrending( ...
            case_name, ig_date, post_date, L8_SI, L8_TIME, MTBS_BNDY, ...
            NLCD_LC, VR, LON, LAT, boundaries, paths, cfg)
% Perform phenology detrending using the VR as the fitting region.
%
% For each pre-fire date
%   1. Use SI_pre and SI_post pairs within the VR
%   2. Fit a linear relationship by land-cover class
%   3. Apply the fitted model to SI_pre to obtain pSI_pre
%   4. Derive conventional and phenology-detrended burn severity indices

    post_idx = find(L8_TIME == post_date, 1, 'first');
    if isempty(post_idx)
        error('Post-fire date (%s) was not found in L8_TIME for case: %s', ...
            datestr(post_date, 'yyyymmdd'), case_name);
    end

    pre_idx = find(L8_TIME <= ig_date);
    if isempty(pre_idx)
        error('No pre-fire images were found before ignition for case: %s', case_name);
    end

    L8_TIME_pre = L8_TIME(pre_idx);
    data_post = L8_SI(:, :, :, post_idx);

    nrow = size(MTBS_BNDY, 1);
    ncol = size(MTBS_BNDY, 2);
    npre = numel(pre_idx);
    nlc = 6;
    nsi = 2;

    L8_dSI  = nan(nrow, ncol, numel(cfg.list_dsi), npre, 'single');
    L8_pdSI = nan(nrow, ncol, numel(cfg.list_pdsi), npre, 'single');

    % Dimensions: pre-date x land-cover-class x SI x [slope, intercept]
    COEF = nan(npre, nlc, nsi, 2, 'single');

    for dt = 1:npre
        data_pre = L8_SI(:, :, :, pre_idx(dt));

        d_maps = nan(nrow, ncol, numel(cfg.list_dsi), 'single');
        pd_maps = nan(nrow, ncol, numel(cfg.list_pdsi), 'single');

        for si = 1:nsi
            X = data_pre(:, :, si);
            Y = data_post(:, :, si);
            pSI = X;

            for lc = 1:nlc
                valid_mask = (VR == 1) & (NLCD_LC == lc);

                X_vr = X(valid_mask);
                Y_vr = Y(valid_mask);

                valid_pairs = ~isnan(X_vr) & ~isnan(Y_vr);
                XY = unique([X_vr(valid_pairs), Y_vr(valid_pairs)], 'rows');

                if size(XY, 1) >= cfg.min_lr_samples
                    ols_coef = polyfit(XY(:, 1), XY(:, 2), 1);
                else
                    ols_coef = [nan, nan];
                end

                pSI(NLCD_LC == lc) = X(NLCD_LC == lc) * ols_coef(1) + ols_coef(2);
                COEF(dt, lc, si, :) = single(ols_coef);
            end

            if si == 1
                d_maps(:, :, 1) = X - Y;
                pd_maps(:, :, 1) = pSI - Y;
            else
                d_maps(:, :, 2) = X - Y;
                pd_maps(:, :, 2) = pSI - Y;

                d_maps(:, :, 3) = (X - Y) ./ sqrt(abs(X));
                pd_maps(:, :, 3) = (pSI - Y) ./ sqrt(abs(pSI));

                d_maps(:, :, 4) = (X - Y) ./ (X + 1.001);
                pd_maps(:, :, 4) = (pSI - Y) ./ (pSI + 1.001);
            end
        end

        L8_dSI(:, :, :, dt) = d_maps;
        L8_pdSI(:, :, :, dt) = pd_maps;

        save_index_map_set(LON, LAT, d_maps, boundaries, case_name, L8_TIME_pre(dt), ...
            cfg.list_dsi, paths.step3, cfg.dsi_clim, 'dsi', cfg);

        save_index_map_set(LON, LAT, pd_maps, boundaries, case_name, L8_TIME_pre(dt), ...
            cfg.list_pdsi, paths.step3, cfg.dsi_clim, 'dsi', cfg);
    end
end

function save_index_map_set(LON, LAT, map_stack, boundaries, case_name, date_value, ...
                            index_names, path_base, clim_range, map_type, cfg)
% Save a set of index maps using the same rendering settings.

    for i = 1:numel(index_names)
        save_index_map(LON, LAT, map_stack(:, :, i), boundaries, case_name, ...
            date_value, index_names{i}, path_base, clim_range, map_type, cfg);
    end
end

function save_fire_map(fig_path, case_name, LON, LAT, MTBS_BNDY, pxid_modis, pxid_viirs, ...
                       af_boundaries, roi_bbox, vr_bbox, cfg)
% Save the step 1 summary map showing
%   - MTBS/WFIGS boundary raster
%   - MODIS/VIIRS active fire detections
%   - Buffered AF boundary
%   - ROI and VR rectangles

    [fig, ax] = create_map_figure(cfg);
    ax.Color = cfg.fire_map_bg_color;

    plot_georaster(ax, LON, LAT, double(MTBS_BNDY > 0));
    colormap(ax, build_fire_map_colormap());
    caxis(ax, [0, 1]);

    modis_legend = scatter(ax, nan, nan, 30, 'o', ...
        'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
    viirs_legend = scatter(ax, nan, nan, 18, 'o', ...
        'MarkerEdgeColor', [0.58, 0, 0.83], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
    buffer_legend = plot(ax, nan, nan, '-', 'Color', [0.88, 0.88, 0.78], 'LineWidth', 1.8);
    roi_legend = plot(ax, nan, nan, 'k-', 'LineWidth', 1.5);
    vr_legend  = plot(ax, nan, nan, 'b-', 'LineWidth', 1.5);

    if ~isempty(pxid_modis)
        scatter(ax, LON(pxid_modis), LAT(pxid_modis), 30, 'o', ...
            'MarkerEdgeColor', [1, 0, 0], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
    end

    if ~isempty(pxid_viirs)
        scatter(ax, LON(pxid_viirs), LAT(pxid_viirs), 18, 'o', ...
            'MarkerEdgeColor', [0.58, 0, 0.83], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
    end

    plot_boundaries(ax, LON, LAT, af_boundaries, [0.88, 0.88, 0.78], 1.8);
    draw_bbox(ax, LON, LAT, roi_bbox, 'k', 1.5);
    draw_bbox(ax, LON, LAT, vr_bbox, 'b', 1.5);

    format_map_axes(ax, LON, LAT, cfg.map_font_size);
    title(ax, case_name, 'Interpreter', 'none', 'FontSize', cfg.step1_title_font_size);

    legend(ax, [modis_legend, viirs_legend, buffer_legend, roi_legend, vr_legend], ...
        {'MODIS AF', 'VIIRS AF', 'Buffered AF', 'ROI', 'VR'}, ...
        'Location', 'best', 'Box', 'on', 'Color', 'w');

    print(fig, fig_path, '-dpng', '-r300');
    close(fig);
end

function save_index_map(LON, LAT, index_map, boundaries, case_name, date_value, ...
                        index_name, path_base, clim_range, map_type, cfg)
% Save a map of a spectral index or burn severity index.
%
% Inputs
%   map_type = 'si'  -> single-date spectral index style
%   map_type = 'dsi' -> differenced index style

    [fig, ax] = create_map_figure(cfg);

    plot_georaster(ax, LON, LAT, index_map);
    colormap(ax, build_index_colormap(map_type));

    cb = colorbar(ax);
    cb.FontSize = cfg.colorbar_font_size;

    caxis(ax, clim_range);
    plot_boundaries(ax, LON, LAT, boundaries, 'r', 2);
    format_map_axes(ax, LON, LAT, cfg.map_font_size);

    file_name = sprintf('%s_%s_%s', case_name, index_name, datestr(date_value, 'yyyymmdd'));
    title(ax, file_name, 'FontSize', cfg.title_font_size, 'Interpreter', 'none');

    outdir = fullfile(path_base, 'Figures', index_name);
    ensure_dir(outdir);

    print(fig, fullfile(outdir, [file_name, '.png']), '-dpng', '-r300');
    close(fig);
end

function [fig, ax] = create_map_figure(cfg)
% Create a common off-screen figure for map rendering.

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', cfg.figure_position);
    ax = axes(fig);
    hold(ax, 'on');
end

function plot_georaster(ax, LON, LAT, raster_data)
% Plot a raster in longitude-latitude coordinates using built-in graphics.

    Z = zeros(size(raster_data));
    alpha_mask = ~isnan(raster_data);

    s = surface(ax, LON, LAT, Z, raster_data, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'texturemap');

    set(s, 'AlphaData', double(alpha_mask), 'FaceAlpha', 'texturemap');
    view(ax, 2);
end

function plot_boundaries(ax, LON, LAT, boundaries, line_color, line_width)
% Plot pixel-based boundaries in longitude-latitude coordinates.

    for k = 1:numel(boundaries)
        B = boundaries{k};
        row_idx = B(:, 1);
        col_idx = B(:, 2);

        lon_b = LON(sub2ind(size(LON), row_idx, col_idx));
        lat_b = LAT(sub2ind(size(LAT), row_idx, col_idx));

        plot(ax, lon_b, lat_b, '-', 'Color', line_color, 'LineWidth', line_width);
    end
end

function draw_bbox(ax, LON, LAT, bbox, line_color, line_width)
% Draw a rectangular bounding box defined in row/column coordinates.

    [lon_box, lat_box] = bbox_to_polygon_coords(LON, LAT, bbox);
    plot(ax, lon_box, lat_box, '-', 'Color', line_color, 'LineWidth', line_width);
end

function [lon_box, lat_box] = bbox_to_polygon_coords(LON, LAT, bbox)
% Convert [row_min, row_max, col_min, col_max] into polygon coordinates.

    rows = [bbox(1), bbox(2), bbox(2), bbox(1), bbox(1)];
    cols = [bbox(3), bbox(3), bbox(4), bbox(4), bbox(3)];

    lon_box = LON(sub2ind(size(LON), rows, cols));
    lat_box = LAT(sub2ind(size(LAT), rows, cols));
end

function format_map_axes(ax, LON, LAT, font_size)
% Apply consistent formatting to lon-lat map axes.

    xlim(ax, [min(LON(:)), max(LON(:))]);
    ylim(ax, [min(LAT(:)), max(LAT(:))]);

    axis(ax, 'tight');
    axis(ax, 'xy');
    box(ax, 'on');
    grid(ax, 'on');

    ax.FontSize = font_size;
    ax.LineWidth = 1;
    ax.Layer = 'top';

    xlabel(ax, 'Longitude', 'FontSize', font_size);
    ylabel(ax, 'Latitude',  'FontSize', font_size);
end

function cmap = build_fire_map_colormap()
% Colormap for the step 1 map.

    cmap = [ ...
        0.24, 0.38, 0.24;
        0.52, 0.69, 0.43];
end

function cmap = build_index_colormap(map_type)
% Colormap for spectral-index and differenced-index maps.

    base = [ ...
        231,  70,  70;
        250, 152, 132;
        248, 196, 180;
        245, 236, 213;
        144, 198, 124;
        103, 174, 110;
         50, 142, 110] / 255;

    cmap_full = interp1(linspace(0, 1, size(base, 1)), ...
                        base, linspace(0, 1, 256));

    switch lower(map_type)
        case 'si'
            cmap = cmap_full(129:end, :);
        case 'dsi'
            cmap = flipud(cmap_full);
        otherwise
            error('Unknown map_type: %s', map_type);
    end
end
