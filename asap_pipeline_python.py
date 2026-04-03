
#!/usr/bin/env python3
"""
ASAP burn severity analysis pipeline in Python.

This script integrates the three ASAP preprocessing steps into a single file:

    Step 1. Definition of the analysis area
    Step 2. Preparation of analysis imagery
    Step 3. Phenology detrending

Each step is saved separately:

    01_ASAP_step1/Dataset/*.nc
    02_ASAP_step2/Dataset/*.nc
    03_ASAP_step3/Dataset/*.nc
"""

from __future__ import annotations

import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.ndimage import binary_dilation, find_objects, label
from skimage.measure import find_contours


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

@dataclass
class Config:
    """User settings and constants."""

    clear_ratio_threshold: float = 0.70
    af_buffer_size: int = 201
    roi_padding_pixels: int = 100
    vr_step_pixels: int = 100
    vr_max_steps: int = 5
    min_lr_samples: int = 10
    reflectance_min_dn: int = 7273
    reflectance_max_dn: int = 43636
    firms_window_days: int = 28

    qa_clear_value: int = 21824
    qa_exclude_value: int = 21952
    reflectance_scale: float = 0.0000275
    reflectance_offset: float = -0.2

    figure_size_inches: Tuple[float, float] = (10.0, 10.0)
    figure_dpi: int = 300
    map_font_size: int = 12
    colorbar_font_size: int = 10
    title_font_size: int = 14
    step1_title_font_size: int = 13
    fire_map_bg_color: Tuple[float, float, float] = (0.18, 0.28, 0.18)
    si_clim: Tuple[float, float] = (0.0, 0.7)
    dsi_clim: Tuple[float, float] = (-0.7, 0.7)

    list_si: Tuple[str, str] = ("NDVI", "NBR")
    list_dsi: Tuple[str, str, str, str] = ("dNDVI", "dNBR", "rdNBR", "RBR")
    list_pdsi: Tuple[str, str, str, str] = ("pdNDVI", "pdNBR", "prdNBR", "pRBR")

    nc_engine: str | None = None
    nc_compression_level: int = 4


# -----------------------------------------------------------------------------
# Path helpers
# -----------------------------------------------------------------------------

def get_project_paths(path_root: Path) -> Dict[str, Path]:
    """Build all project paths from the project root."""
    step1 = path_root / "01_ASAP_step1"
    step2 = path_root / "02_ASAP_step2"
    step3 = path_root / "03_ASAP_step3"

    return {
        "root": path_root,
        "data": path_root / "00_Data",
        "step1": step1,
        "step2": step2,
        "step3": step3,
        "step1_dataset": step1 / "Dataset",
        "step1_figures": step1 / "Figures",
        "step2_dataset": step2 / "Dataset",
        "step2_figures": step2 / "Figures",
        "step3_dataset": step3 / "Dataset",
        "step3_figures": step3 / "Figures",
    }


def ensure_dir(folder_path: Path) -> None:
    """Create a folder if it does not already exist."""
    folder_path.mkdir(parents=True, exist_ok=True)


def create_output_dirs(paths: Dict[str, Path]) -> None:
    """Create all step output directories."""
    keys = [
        "step1_dataset", "step1_figures",
        "step2_dataset", "step2_figures",
        "step3_dataset", "step3_figures",
    ]
    for key in keys:
        ensure_dir(paths[key])


# -----------------------------------------------------------------------------
# Raster I/O and preprocessing
# -----------------------------------------------------------------------------

def build_lon_lat_grids(transform: rasterio.Affine, height: int, width: int) -> Tuple[np.ndarray, np.ndarray]:
    """Build 2D longitude and latitude grids from a raster transform."""
    cols = np.arange(width)
    rows = np.arange(height)

    xs = np.array(rasterio.transform.xy(transform, np.zeros_like(cols), cols, offset="center")[0], dtype=np.float64)
    ys = np.array(rasterio.transform.xy(transform, rows, np.zeros_like(rows), offset="center")[1], dtype=np.float64)

    lon, lat = np.meshgrid(xs, ys)
    return lon, lat


def read_single_band_tif(tif_path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a single-band GeoTIFF and return data, lon grid, and lat grid."""
    with rasterio.open(tif_path) as src:
        data = src.read(1)
        lon, lat = build_lon_lat_grids(src.transform, src.height, src.width)
    return data, lon, lat


def read_mtbs_boundary(case_name: str, path_data: Path):
    folder = path_data / "MTBS_WFIGS"

    # 1) exact match with .tif
    tif_path = folder / f"{case_name}.tif"
    if tif_path.exists():
        return read_single_band_tif(tif_path)

    # 2) fallback: search similar tif names
    matches = sorted(folder.glob(f"{case_name}*.tif"))
    if not matches:
        raise FileNotFoundError(
            f"MTBS/WFIGS tif not found for case '{case_name}' in {folder}"
        )

    return read_single_band_tif(matches[0])


def read_nlcd_reclassified(case_name: str, path_data: Path) -> np.ndarray:
    """
    Read the NLCD raster for one wildfire case and reclassify it into the
    six vegetation classes used in ASAP.

    Reclassified classes
        0: Excluded / non-vegetated
        1: Deciduous forest
        2: Evergreen forest
        3: Mixed forest
        4: Shrubland
        5: Herbaceous
        6: Planted/cultivated
    """
    folder_nlcd = path_data / "NLCD"
    matches = sorted(folder_nlcd.glob(f"{case_name}_*.tif"))

    if not matches:
        raise FileNotFoundError(f"NLCD file not found for case: {case_name}")

    with rasterio.open(matches[0]) as src:
        nlcd = src.read(1)

    nlcd = nlcd.astype(np.int16, copy=False)
    nlcd[(nlcd < 41) | (nlcd > 82)] = 0
    nlcd[nlcd == 41] = 1
    nlcd[nlcd == 42] = 2
    nlcd[nlcd == 43] = 3
    nlcd[(nlcd >= 51) & (nlcd <= 52)] = 4
    nlcd[(nlcd >= 71) & (nlcd <= 74)] = 5
    nlcd[(nlcd >= 81) & (nlcd <= 82)] = 6

    return nlcd.astype(np.uint8)


def rasterize_firms_detections(
    firms_table: pd.DataFrame,
    ig_date: pd.Timestamp,
    lon: np.ndarray,
    lat: np.ndarray,
    nlcd_lc: np.ndarray,
    cfg: Config,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert FIRMS point detections into a raster aligned with the MTBS/NLCD grid.

    Processing steps
        1. Keep detections inside the temporal range
        2. Keep detections inside the target image extent
        3. Assign each detection to the nearest raster cell
        4. Remove detections over excluded land-cover classes
    """
    ig_date = pd.Timestamp(ig_date).normalize()
    time_mask = (firms_table["acq_date"] >= ig_date) & (
        firms_table["acq_date"] <= ig_date + pd.Timedelta(days=cfg.firms_window_days)
    )
    firms_sub = firms_table.loc[time_mask].copy()

    lat_min, lat_max = np.nanmin(lat), np.nanmax(lat)
    lon_min, lon_max = np.nanmin(lon), np.nanmax(lon)

    firms_sub = firms_sub.loc[
        (firms_sub["latitude"] >= lat_min)
        & (firms_sub["latitude"] <= lat_max)
        & (firms_sub["longitude"] >= lon_min)
        & (firms_sub["longitude"] <= lon_max)
    ]

    af_map = np.zeros(lon.shape, dtype=np.uint8)

    if firms_sub.empty:
        return af_map, np.array([], dtype=np.int64)

    lon_axis = lon[0, :]
    lat_axis = lat[:, 0]

    for _, row in firms_sub.iterrows():
        col_idx = int(np.argmin(np.abs(lon_axis - row["longitude"])))
        row_idx = int(np.argmin(np.abs(lat_axis - row["latitude"])))
        af_map[row_idx, col_idx] = 1

    af_map[nlcd_lc == 0] = 0
    pxid = np.flatnonzero(af_map == 1)

    return af_map, pxid


def connected_component_bboxes(mask: np.ndarray) -> List[Tuple[int, int, int, int]]:
    """
    Convert connected components in a binary mask into bounding boxes:
    (row_min, row_max, col_min, col_max), using inclusive indexing.
    """
    labeled, _ = label(mask.astype(bool))
    slices = find_objects(labeled)

    bboxes: List[Tuple[int, int, int, int]] = []
    for slc in slices:
        if slc is None:
            continue
        row_slice, col_slice = slc
        bboxes.append(
            (
                int(row_slice.start),
                int(row_slice.stop - 1),
                int(col_slice.start),
                int(col_slice.stop - 1),
            )
        )
    return bboxes


def select_center_roi_bbox(roi_bboxes: Sequence[Tuple[int, int, int, int]], img_shape: Tuple[int, int]) -> int:
    """
    Select the ROI candidate whose bounding box contains the image center.

    Returns the index of the selected bbox. If no candidate contains the
    center, return -1 and use the fallback ROI.
    """
    if not roi_bboxes:
        return -1

    center_row = (img_shape[0] - 1) / 2.0
    center_col = (img_shape[1] - 1) / 2.0

    for i, (row_min, row_max, col_min, col_max) in enumerate(roi_bboxes):
        if row_min <= center_row <= row_max and col_min <= center_col <= col_max:
            return i

    return -1


def compute_mask_bbox(mask: np.ndarray, padding: int) -> Tuple[int, int, int, int] | None:
    """
    Convert a binary mask into a padded bounding box:
    (row_min, row_max, col_min, col_max), inclusive.
    """
    rows, cols = np.where(mask)
    if rows.size == 0:
        return None

    row_min = max(int(rows.min()) - padding, 0)
    row_max = min(int(rows.max()) + padding, mask.shape[0] - 1)
    col_min = max(int(cols.min()) - padding, 0)
    col_max = min(int(cols.max()) + padding, mask.shape[1] - 1)

    return row_min, row_max, col_min, col_max


def mask_to_bbox(mask: np.ndarray) -> Tuple[int, int, int, int]:
    """Convert a binary mask into its tight bounding box."""
    bbox = compute_mask_bbox(mask.astype(bool), padding=0)
    if bbox is None:
        raise ValueError("Bounding box cannot be generated because the mask is empty.")
    return bbox


def bbox_to_slices(bbox: Tuple[int, int, int, int]) -> Tuple[slice, slice]:
    """Convert an inclusive bbox into NumPy slices."""
    row_min, row_max, col_min, col_max = bbox
    return slice(row_min, row_max + 1), slice(col_min, col_max + 1)


def define_roi(
    mtbs_bndy: np.ndarray,
    roi_bboxes: Sequence[Tuple[int, int, int, int]],
    target_idx: int,
    roi_padding_pixels: int,
) -> Tuple[np.ndarray, np.ndarray, Tuple[int, int, int, int]]:
    """
    Define the rectangular ROI.

    If a target AF-based ROI exists, use that bounding box.
    Otherwise, build a fallback ROI from the MTBS burn perimeter with padding.
    """
    roi = np.zeros(mtbs_bndy.shape, dtype=np.uint8)
    roi_af = np.zeros(mtbs_bndy.shape, dtype=np.uint8)

    if target_idx < 0:
        target_bbox = compute_mask_bbox(mtbs_bndy == 1, roi_padding_pixels)
        if target_bbox is None:
            raise ValueError("Fallback ROI cannot be defined because MTBS_BNDY is empty.")
    else:
        target_bbox = roi_bboxes[target_idx]

    roi_rows, roi_cols = bbox_to_slices(target_bbox)
    roi[roi_rows, roi_cols] = 1

    for bbox in roi_bboxes:
        rows, cols = bbox_to_slices(bbox)
        roi_af[rows, cols] = 1

    return roi, roi_af, target_bbox


def define_vr(
    roi: np.ndarray,
    roi_af: np.ndarray,
    nlcd_lc: np.ndarray,
    target_bbox: Tuple[int, int, int, int],
    vr_step_pixels: int,
    vr_max_steps: int,
) -> np.ndarray:
    """
    Define the vicinity region (VR) by progressively expanding a rectangle
    around the target ROI until all vegetation classes in the ROI are also
    represented in the VR.

    The VR excludes all AF-based ROI candidates and the target ROI itself.
    """
    vr = np.zeros(roi.shape, dtype=np.uint8)

    roi_classes = np.unique(nlcd_lc[roi == 1])
    roi_classes = roi_classes[roi_classes != 0]

    for j in range(1, vr_max_steps + 1):
        row_min = max(target_bbox[0] - vr_step_pixels * j, 0)
        row_max = min(target_bbox[1] + vr_step_pixels * j, roi.shape[0] - 1)
        col_min = max(target_bbox[2] - vr_step_pixels * j, 0)
        col_max = min(target_bbox[3] + vr_step_pixels * j, roi.shape[1] - 1)

        vr.fill(0)
        vr[row_min:row_max + 1, col_min:col_max + 1] = 1

        vr[roi_af == 1] = 0
        vr[roi == 1] = 0

        vr_classes = np.unique(nlcd_lc[vr == 1])
        vr_classes = vr_classes[vr_classes != 0]

        if np.all(np.isin(roi_classes, vr_classes)):
            break

    return vr


def crop_base_layers(
    lat: np.ndarray,
    lon: np.ndarray,
    mtbs_bndy: np.ndarray,
    nlcd_lc: np.ndarray,
    roi: np.ndarray,
    vr: np.ndarray,
    crop_bbox: Tuple[int, int, int, int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Crop all base layers using the same bounding box."""
    rows, cols = bbox_to_slices(crop_bbox)

    return (
        lat[rows, cols],
        lon[rows, cols],
        mtbs_bndy[rows, cols],
        nlcd_lc[rows, cols],
        roi[rows, cols],
        vr[rows, cols],
    )


# -----------------------------------------------------------------------------
# Landsat processing
# -----------------------------------------------------------------------------

def extract_date_from_name(file_name: str) -> pd.Timestamp:
    """Extract YYYYMMDD from the last 8 characters before the file extension."""
    return pd.to_datetime(file_name[-12:-4], format="%Y%m%d")


def compute_clear_ratio(qa_band: np.ndarray, cfg: Config) -> float:
    """Compute the clear ratio used in the original pipeline."""
    numerator = np.sum(qa_band == cfg.qa_clear_value)
    denominator = np.sum((qa_band > 1) & (qa_band != cfg.qa_exclude_value))

    if denominator == 0:
        return 0.0
    return float(numerator) / float(denominator)


def apply_landsat_qc(
    reflectance_bands: np.ndarray,
    qa_band: np.ndarray,
    nlcd_lc: np.ndarray,
    cfg: Config,
) -> np.ndarray:
    """Apply Landsat QA filtering, land-cover masking, and DN-to-reflectance conversion."""
    out = reflectance_bands.astype(np.float32, copy=True)

    qa_mask = qa_band != cfg.qa_clear_value
    lc_mask = nlcd_lc == 0

    out[qa_mask, :] = np.nan
    out[lc_mask, :] = np.nan

    invalid_dn = (out < cfg.reflectance_min_dn) | (out > cfg.reflectance_max_dn)
    out[invalid_dn] = np.nan

    out = out * cfg.reflectance_scale + cfg.reflectance_offset
    return out


def calculate_spectral_indices(reflectance_bands: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate NDVI and NBR from Landsat surface reflectance bands.

    Landsat band order used here
        band 3: red
        band 4: NIR
        band 6: SWIR2
    """
    red = reflectance_bands[:, :, 2]
    nir = reflectance_bands[:, :, 3]
    swir2 = reflectance_bands[:, :, 5]

    with np.errstate(invalid="ignore", divide="ignore"):
        ndvi = (nir - red) / (nir + red)
        nbr = (nir - swir2) / (nir + swir2)

    return ndvi.astype(np.float32), nbr.astype(np.float32)


def prepare_landsat_spectral_indices(
    case_name: str,
    pre_date: pd.Timestamp,
    post_date: pd.Timestamp,
    paths: Dict[str, Path],
    crop_bbox: Tuple[int, int, int, int],
    lon: np.ndarray,
    lat: np.ndarray,
    nlcd_lc: np.ndarray,
    burn_boundaries: Sequence[np.ndarray],
    cfg: Config,
) -> Tuple[np.ndarray, pd.DatetimeIndex]:
    """
    Read Landsat scenes for one case, apply QA and land-cover filtering,
    compute NDVI and NBR, and save a figure for every retained scene.

    Retention rule
        - Always keep the metadata-based pre- and post-fire scenes
        - For all other dates, keep only scenes whose clear ratio exceeds the
          specified threshold
    """
    list_l8 = sorted((paths["data"] / "Landsat").glob(f"{case_name}*.tif"), key=lambda p: extract_date_from_name(p.name))

    if not list_l8:
        raise FileNotFoundError(f"No Landsat files found for case: {case_name}")

    si_list: List[np.ndarray] = []
    time_list: List[pd.Timestamp] = []

    rows, cols = bbox_to_slices(crop_bbox)
    pre_date = pd.Timestamp(pre_date).normalize()
    post_date = pd.Timestamp(post_date).normalize()

    for tif_path in list_l8:
        img_date = extract_date_from_name(tif_path.name).normalize()

        with rasterio.open(tif_path) as src:
            l8_raw = src.read().transpose(1, 2, 0).astype(np.float32)

        l8_raw = l8_raw[rows, cols, :]

        reflectance_bands = l8_raw[:, :, :6]
        qa_band = l8_raw[:, :, -1]

        keep_scene = (img_date == post_date) or (img_date == pre_date)
        if not keep_scene:
            clear_ratio = compute_clear_ratio(qa_band, cfg)
            if clear_ratio < cfg.clear_ratio_threshold:
                continue

        reflectance_bands = apply_landsat_qc(reflectance_bands, qa_band, nlcd_lc, cfg)
        ndvi, nbr = calculate_spectral_indices(reflectance_bands)

        si_list.append(np.stack([ndvi, nbr], axis=2))
        time_list.append(img_date)

        save_index_map(
            lon, lat, ndvi, burn_boundaries, case_name, img_date,
            "NDVI", paths["step2"], cfg.si_clim, "si", cfg
        )
        save_index_map(
            lon, lat, nbr, burn_boundaries, case_name, img_date,
            "NBR", paths["step2"], cfg.si_clim, "si", cfg
        )

    if not si_list:
        raise ValueError(f"No Landsat scenes were retained for case: {case_name}")

    l8_si = np.stack(si_list, axis=3).astype(np.float32)
    l8_time = pd.DatetimeIndex(time_list)

    return l8_si, l8_time


# -----------------------------------------------------------------------------
# Phenology detrending
# -----------------------------------------------------------------------------

def run_phenology_detrending(
    case_name: str,
    ig_date: pd.Timestamp,
    post_date: pd.Timestamp,
    l8_si: np.ndarray,
    l8_time: pd.DatetimeIndex,
    mtbs_bndy: np.ndarray,
    nlcd_lc: np.ndarray,
    vr: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    boundaries: Sequence[np.ndarray],
    paths: Dict[str, Path],
    cfg: Config,
) -> Tuple[np.ndarray, np.ndarray, pd.DatetimeIndex, np.ndarray]:
    """
    Perform phenology detrending using the VR as the fitting region.

    For each pre-fire date
        1. Use SI_pre and SI_post pairs within the VR
        2. Fit a linear relationship by land-cover class
        3. Apply the fitted model to SI_pre to obtain pSI_pre
        4. Derive conventional and phenology-detrended burn severity indices
    """
    post_date = pd.Timestamp(post_date).normalize()
    ig_date = pd.Timestamp(ig_date).normalize()

    post_matches = np.where(l8_time.normalize() == post_date)[0]
    if post_matches.size == 0:
        raise ValueError(f'Post-fire date ({post_date.strftime("%Y%m%d")}) was not found in L8_TIME for case: {case_name}')
    post_idx = int(post_matches[0])

    pre_idx = np.where(l8_time.normalize() <= ig_date)[0]
    if pre_idx.size == 0:
        raise ValueError(f"No pre-fire images were found before ignition for case: {case_name}")

    l8_time_pre = pd.DatetimeIndex(l8_time[pre_idx])
    data_post = l8_si[:, :, :, post_idx]

    nrow, ncol = mtbs_bndy.shape
    npre = pre_idx.size
    nlc = 6
    nsi = 2

    l8_dsi = np.full((nrow, ncol, len(cfg.list_dsi), npre), np.nan, dtype=np.float32)
    l8_pdsi = np.full((nrow, ncol, len(cfg.list_pdsi), npre), np.nan, dtype=np.float32)
    coef = np.full((npre, nlc, nsi, 2), np.nan, dtype=np.float32)

    for dt, src_idx in enumerate(pre_idx):
        data_pre = l8_si[:, :, :, src_idx]

        d_maps = np.full((nrow, ncol, len(cfg.list_dsi)), np.nan, dtype=np.float32)
        pd_maps = np.full((nrow, ncol, len(cfg.list_pdsi)), np.nan, dtype=np.float32)

        for si in range(nsi):
            x = data_pre[:, :, si]
            y = data_post[:, :, si]
            psi = x.copy()

            for lc in range(1, nlc + 1):
                valid_mask = (vr == 1) & (nlcd_lc == lc)

                x_vr = x[valid_mask]
                y_vr = y[valid_mask]

                valid_pairs = ~np.isnan(x_vr) & ~np.isnan(y_vr)
                xy = np.column_stack([x_vr[valid_pairs], y_vr[valid_pairs]])

                if xy.size == 0:
                    fit_coef = np.array([np.nan, np.nan], dtype=np.float32)
                else:
                    xy = np.unique(xy, axis=0)
                    if xy.shape[0] >= cfg.min_lr_samples:
                        fit_coef = np.polyfit(xy[:, 0], xy[:, 1], deg=1).astype(np.float32)
                    else:
                        fit_coef = np.array([np.nan, np.nan], dtype=np.float32)

                psi[nlcd_lc == lc] = x[nlcd_lc == lc] * fit_coef[0] + fit_coef[1]
                coef[dt, lc - 1, si, :] = fit_coef

            if si == 0:
                d_maps[:, :, 0] = x - y
                pd_maps[:, :, 0] = psi - y
            else:
                with np.errstate(invalid="ignore", divide="ignore"):
                    d_maps[:, :, 1] = x - y
                    pd_maps[:, :, 1] = psi - y

                    d_maps[:, :, 2] = (x - y) / np.sqrt(np.abs(x))
                    pd_maps[:, :, 2] = (psi - y) / np.sqrt(np.abs(psi))

                    d_maps[:, :, 3] = (x - y) / (x + 1.001)
                    pd_maps[:, :, 3] = (psi - y) / (psi + 1.001)

        l8_dsi[:, :, :, dt] = d_maps
        l8_pdsi[:, :, :, dt] = pd_maps

        save_index_map_set(
            lon, lat, d_maps, boundaries, case_name, l8_time_pre[dt],
            cfg.list_dsi, paths["step3"], cfg.dsi_clim, "dsi", cfg
        )
        save_index_map_set(
            lon, lat, pd_maps, boundaries, case_name, l8_time_pre[dt],
            cfg.list_pdsi, paths["step3"], cfg.dsi_clim, "dsi", cfg
        )

    return l8_dsi, l8_pdsi, l8_time_pre, coef


# -----------------------------------------------------------------------------
# Plotting helpers
# -----------------------------------------------------------------------------

def mask_to_contours(mask: np.ndarray) -> List[np.ndarray]:
    """
    Convert a binary mask into contour coordinates in pixel space.

    Returns a list of arrays with shape (N, 2), where columns are
    [row, col] in floating-point pixel coordinates.
    """
    contours = find_contours(mask.astype(float), level=0.5)
    return [np.asarray(c, dtype=np.float64) for c in contours if c.shape[0] >= 2]


def contour_pixels_to_lonlat(
    contour_rc: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert floating-point pixel contour coordinates into lon/lat coordinates."""
    row_coords = contour_rc[:, 0]
    col_coords = contour_rc[:, 1]

    lon_axis = lon[0, :]
    lat_axis = lat[:, 0]

    lon_line = np.interp(col_coords, np.arange(lon.shape[1]), lon_axis)
    lat_line = np.interp(row_coords, np.arange(lat.shape[0]), lat_axis)

    return lon_line, lat_line


def bbox_to_polygon_coords(
    lon: np.ndarray,
    lat: np.ndarray,
    bbox: Tuple[int, int, int, int],
) -> Tuple[np.ndarray, np.ndarray]:
    """Convert an inclusive bbox into polygon coordinates."""
    row_min, row_max, col_min, col_max = bbox

    rows = np.array([row_min, row_max, row_max, row_min, row_min], dtype=int)
    cols = np.array([col_min, col_min, col_max, col_max, col_min], dtype=int)

    lon_box = lon[rows, cols]
    lat_box = lat[rows, cols]

    return lon_box, lat_box


def create_map_figure(cfg: Config) -> Tuple[plt.Figure, plt.Axes]:
    """Create a common off-screen figure for map rendering."""
    fig, ax = plt.subplots(figsize=cfg.figure_size_inches, dpi=cfg.figure_dpi)
    return fig, ax


def build_fire_map_colormap() -> ListedColormap:
    """Colormap for the step 1 fire map."""
    colors = np.array([
        [0.24, 0.38, 0.24],
        [0.52, 0.69, 0.43],
    ])
    return ListedColormap(colors)


def build_index_colormap(map_type: str) -> ListedColormap:
    """Colormap for spectral-index and differenced-index maps."""
    base = np.array([
        [231,  70,  70],
        [250, 152, 132],
        [248, 196, 180],
        [245, 236, 213],
        [144, 198, 124],
        [103, 174, 110],
        [ 50, 142, 110],
    ], dtype=np.float32) / 255.0

    cmap_full = LinearSegmentedColormap.from_list("asap_index", base, N=256)
    colors = cmap_full(np.linspace(0, 1, 256))

    if map_type.lower() == "si":
        return ListedColormap(colors[128:, :])
    if map_type.lower() == "dsi":
        return ListedColormap(colors[::-1, :])

    raise ValueError(f"Unknown map_type: {map_type}")


def plot_georaster(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray, raster_data: np.ndarray, cmap) -> None:
    """Plot a raster in lon/lat coordinates."""
    masked = np.ma.masked_invalid(raster_data)
    ax.pcolormesh(lon, lat, masked, shading="auto", cmap=cmap)


def plot_contours(
    ax: plt.Axes,
    lon: np.ndarray,
    lat: np.ndarray,
    contours: Sequence[np.ndarray],
    line_color,
    line_width: float,
) -> None:
    """Plot contours defined in pixel coordinates on lon/lat axes."""
    for contour in contours:
        lon_line, lat_line = contour_pixels_to_lonlat(contour, lon, lat)
        ax.plot(lon_line, lat_line, "-", color=line_color, linewidth=line_width)


def draw_bbox(
    ax: plt.Axes,
    lon: np.ndarray,
    lat: np.ndarray,
    bbox: Tuple[int, int, int, int],
    line_color,
    line_width: float,
) -> None:
    """Draw a rectangular bounding box defined in row/column coordinates."""
    lon_box, lat_box = bbox_to_polygon_coords(lon, lat, bbox)
    ax.plot(lon_box, lat_box, "-", color=line_color, linewidth=line_width)


def format_map_axes(ax: plt.Axes, lon: np.ndarray, lat: np.ndarray, font_size: int) -> None:
    """Apply consistent formatting to map axes."""
    ax.set_xlim(float(np.nanmin(lon)), float(np.nanmax(lon)))
    ax.set_ylim(float(np.nanmin(lat)), float(np.nanmax(lat)))
    ax.grid(True)
    ax.set_xlabel("Longitude", fontsize=font_size)
    ax.set_ylabel("Latitude", fontsize=font_size)
    ax.tick_params(labelsize=font_size)
    ax.set_aspect("equal", adjustable="box")


def save_fire_map(
    fig_path: Path,
    case_name: str,
    lon: np.ndarray,
    lat: np.ndarray,
    mtbs_bndy: np.ndarray,
    pxid_modis: np.ndarray,
    pxid_viirs: np.ndarray,
    af_contours: Sequence[np.ndarray],
    roi_bbox: Tuple[int, int, int, int],
    vr_bbox: Tuple[int, int, int, int],
    cfg: Config,
) -> None:
    """
    Save the step 1 summary map showing
        - MTBS/WFIGS boundary raster
        - MODIS/VIIRS active fire detections
        - Buffered AF boundary
        - ROI and VR rectangles
    """
    fig, ax = create_map_figure(cfg)
    ax.set_facecolor(cfg.fire_map_bg_color)

    ax.pcolormesh(
        lon, lat, np.ma.masked_invalid((mtbs_bndy > 0).astype(float)),
        shading="auto",
        cmap=build_fire_map_colormap(),
        vmin=0,
        vmax=1,
    )

    modis_handle = ax.scatter([], [], s=30, facecolors="none", edgecolors=(1, 0, 0), linewidths=1.5)
    viirs_handle = ax.scatter([], [], s=18, facecolors="none", edgecolors=(0.58, 0, 0.83), linewidths=1.5)
    buffer_handle, = ax.plot([], [], "-", color=(0.88, 0.88, 0.78), linewidth=1.8)
    roi_handle, = ax.plot([], [], "k-", linewidth=1.5)
    vr_handle, = ax.plot([], [], "b-", linewidth=1.5)

    if pxid_modis.size > 0:
        ax.scatter(
            lon.ravel()[pxid_modis], lat.ravel()[pxid_modis],
            s=30, facecolors="none", edgecolors=(1, 0, 0), linewidths=1.5
        )

    if pxid_viirs.size > 0:
        ax.scatter(
            lon.ravel()[pxid_viirs], lat.ravel()[pxid_viirs],
            s=18, facecolors="none", edgecolors=(0.58, 0, 0.83), linewidths=1.5
        )

    plot_contours(ax, lon, lat, af_contours, (0.88, 0.88, 0.78), 1.8)
    draw_bbox(ax, lon, lat, roi_bbox, "k", 1.5)
    draw_bbox(ax, lon, lat, vr_bbox, "b", 1.5)

    format_map_axes(ax, lon, lat, cfg.map_font_size)
    ax.set_title(case_name, fontsize=cfg.step1_title_font_size)

    ax.legend(
        [modis_handle, viirs_handle, buffer_handle, roi_handle, vr_handle],
        ["MODIS AF", "VIIRS AF", "Buffered AF", "ROI", "VR"],
        loc="best",
        frameon=True,
        facecolor="white",
    )

    fig.tight_layout()
    fig.savefig(fig_path, dpi=cfg.figure_dpi, bbox_inches="tight")
    plt.close(fig)


def save_index_map(
    lon: np.ndarray,
    lat: np.ndarray,
    index_map: np.ndarray,
    boundaries: Sequence[np.ndarray],
    case_name: str,
    date_value: pd.Timestamp,
    index_name: str,
    path_base: Path,
    clim_range: Tuple[float, float],
    map_type: str,
    cfg: Config,
) -> None:
    """Save a map of a spectral index or burn severity index."""
    fig, ax = create_map_figure(cfg)

    mesh = ax.pcolormesh(
        lon, lat, np.ma.masked_invalid(index_map),
        shading="auto",
        cmap=build_index_colormap(map_type),
        vmin=clim_range[0],
        vmax=clim_range[1],
    )
    cb = fig.colorbar(mesh, ax=ax)
    cb.ax.tick_params(labelsize=cfg.colorbar_font_size)

    plot_contours(ax, lon, lat, boundaries, "r", 2.0)
    format_map_axes(ax, lon, lat, cfg.map_font_size)

    file_name = f"{case_name}_{index_name}_{pd.Timestamp(date_value).strftime('%Y%m%d')}"
    ax.set_title(file_name, fontsize=cfg.title_font_size)

    outdir = path_base / "Figures" / index_name
    ensure_dir(outdir)

    fig.tight_layout()
    fig.savefig(outdir / f"{file_name}.png", dpi=cfg.figure_dpi, bbox_inches="tight")
    plt.close(fig)


def save_index_map_set(
    lon: np.ndarray,
    lat: np.ndarray,
    map_stack: np.ndarray,
    boundaries: Sequence[np.ndarray],
    case_name: str,
    date_value: pd.Timestamp,
    index_names: Sequence[str],
    path_base: Path,
    clim_range: Tuple[float, float],
    map_type: str,
    cfg: Config,
) -> None:
    """Save a set of index maps using the same rendering settings."""
    for i, index_name in enumerate(index_names):
        save_index_map(
            lon, lat, map_stack[:, :, i], boundaries, case_name,
            date_value, index_name, path_base, clim_range, map_type, cfg
        )


# -----------------------------------------------------------------------------
# NetCDF output helpers
# -----------------------------------------------------------------------------

def build_netcdf_encoding(ds: xr.Dataset, cfg: Config) -> Dict[str, Dict[str, object]]:
    """Build variable-wise NetCDF encoding for engines that support compression."""
    encoding: Dict[str, Dict[str, object]] = {}
    for var_name in ds.data_vars:
        encoding[var_name] = {
            "zlib": True,
            "complevel": cfg.nc_compression_level,
        }
    return encoding


def resolve_netcdf_engine(preferred_engine: str | None) -> Tuple[str | None, bool]:
    """
    Resolve a NetCDF engine.

    Returns
        engine, supports_compression
    """
    if preferred_engine is not None:
        return preferred_engine, preferred_engine in {"netcdf4", "h5netcdf"}

    try:
        import netCDF4  # noqa: F401
        return "netcdf4", True
    except Exception:
        pass

    try:
        import h5netcdf  # noqa: F401
        return "h5netcdf", True
    except Exception:
        pass

    return None, False


def save_netcdf_dataset(ds: xr.Dataset, out_path: Path, cfg: Config) -> None:
    """Save an xarray Dataset as NetCDF."""
    ensure_dir(out_path.parent)
    engine, supports_compression = resolve_netcdf_engine(cfg.nc_engine)

    if supports_compression:
        encoding = build_netcdf_encoding(ds, cfg)
        if engine is None:
            ds.to_netcdf(out_path, encoding=encoding)
        else:
            ds.to_netcdf(out_path, engine=engine, encoding=encoding)
    else:
        if engine is None:
            ds.to_netcdf(out_path)
        else:
            ds.to_netcdf(out_path, engine=engine)


def make_common_coords(lon: np.ndarray, lat: np.ndarray) -> Dict[str, object]:
    """Build common xarray coordinates."""
    return {
        "y": np.arange(lon.shape[0], dtype=np.int32),
        "x": np.arange(lon.shape[1], dtype=np.int32),
        "lon": (("y", "x"), lon.astype(np.float32)),
        "lat": (("y", "x"), lat.astype(np.float32)),
    }


def save_step1_dataset(
    out_path: Path,
    case_name: str,
    mtbs_bndy: np.ndarray,
    nlcd_lc: np.ndarray,
    roi: np.ndarray,
    vr: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    cfg: Config,
) -> None:
    """Save the step 1 dataset as NetCDF."""
    ds = xr.Dataset(
        data_vars={
            "MTBS_BNDY": (("y", "x"), mtbs_bndy.astype(np.uint8)),
            "NLCD_LC": (("y", "x"), nlcd_lc.astype(np.uint8)),
            "ROI": (("y", "x"), roi.astype(np.uint8)),
            "VR": (("y", "x"), vr.astype(np.uint8)),
        },
        coords=make_common_coords(lon, lat),
        attrs={"case_name": case_name, "storage_format": "NetCDF"},
    )
    save_netcdf_dataset(ds, out_path, cfg)


def save_step2_dataset(
    out_path: Path,
    case_name: str,
    l8_si: np.ndarray,
    l8_time: pd.DatetimeIndex,
    mtbs_bndy: np.ndarray,
    nlcd_lc: np.ndarray,
    roi: np.ndarray,
    vr: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    cfg: Config,
) -> None:
    """Save the step 2 dataset as NetCDF."""
    ds = xr.Dataset(
        data_vars={
            "NDVI": (("y", "x", "time"), l8_si[:, :, 0, :].astype(np.float32)),
            "NBR": (("y", "x", "time"), l8_si[:, :, 1, :].astype(np.float32)),
            "MTBS_BNDY": (("y", "x"), mtbs_bndy.astype(np.uint8)),
            "NLCD_LC": (("y", "x"), nlcd_lc.astype(np.uint8)),
            "ROI": (("y", "x"), roi.astype(np.uint8)),
            "VR": (("y", "x"), vr.astype(np.uint8)),
        },
        coords={
            **make_common_coords(lon, lat),
            "time": pd.DatetimeIndex(l8_time),
        },
        attrs={"case_name": case_name, "storage_format": "NetCDF"},
    )
    save_netcdf_dataset(ds, out_path, cfg)


def save_step3_dataset(
    out_path: Path,
    case_name: str,
    l8_dsi: np.ndarray,
    l8_pdsi: np.ndarray,
    l8_time_pre: pd.DatetimeIndex,
    coef: np.ndarray,
    mtbs_bndy: np.ndarray,
    nlcd_lc: np.ndarray,
    roi: np.ndarray,
    vr: np.ndarray,
    lon: np.ndarray,
    lat: np.ndarray,
    cfg: Config,
) -> None:
    """Save the step 3 dataset as NetCDF."""
    ds = xr.Dataset(
        data_vars={
            "dSI": (("y", "x", "index", "pre_time"), l8_dsi.astype(np.float32)),
            "pdSI": (("y", "x", "index", "pre_time"), l8_pdsi.astype(np.float32)),
            "COEF": (("pre_time", "lc", "si", "param"), coef.astype(np.float32)),
            "MTBS_BNDY": (("y", "x"), mtbs_bndy.astype(np.uint8)),
            "NLCD_LC": (("y", "x"), nlcd_lc.astype(np.uint8)),
            "ROI": (("y", "x"), roi.astype(np.uint8)),
            "VR": (("y", "x"), vr.astype(np.uint8)),
        },
        coords={
            **make_common_coords(lon, lat),
            "index": list(cfg.list_dsi),
            "pre_time": pd.DatetimeIndex(l8_time_pre),
            "lc": np.arange(1, 7, dtype=np.int16),
            "si": list(cfg.list_si),
            "param": ["slope", "intercept"],
        },
        attrs={"case_name": case_name, "storage_format": "NetCDF"},
    )
    save_netcdf_dataset(ds, out_path, cfg)


# -----------------------------------------------------------------------------
# Main pipeline
# -----------------------------------------------------------------------------

def main() -> None:
    """Run the ASAP pipeline."""
    cfg = Config()

    path_root = Path(__file__).resolve().parent
    paths = get_project_paths(path_root)
    create_output_dirs(paths)

    metadata = pd.read_excel(paths["data"] / "metadata.xlsx")
    metadata["Ig_Date"] = pd.to_datetime(metadata["Ig_Date"])
    metadata["Post_Date"] = pd.to_datetime(metadata["Post_Date"])
    metadata["Pre_Date(const)"] = pd.to_datetime(metadata["Pre_Date(const)"])

    firms_modis = pd.read_csv(paths["data"] / "FIRMS" / "fire_M-C61.csv", parse_dates=["acq_date"])
    firms_viirs = pd.read_csv(paths["data"] / "FIRMS" / "fire_V-C2.csv", parse_dates=["acq_date"])

    n_case = len(metadata)

    for ii, meta in metadata.iterrows():
        case_name = str(meta["Incid_Name"])
        ig_date = pd.Timestamp(meta["Ig_Date"])
        post_date = pd.Timestamp(meta["Post_Date"])
        pre_date = pd.Timestamp(meta["Pre_Date(const)"])

        print(f"\n[{ii + 1}/{n_case}] Processing case: {case_name}")

        try:
            # -----------------------------------------------------------------
            # STEP 1. Definition of the analysis area
            # -----------------------------------------------------------------
            mtbs_bndy, lon, lat = read_mtbs_boundary(case_name, paths["data"])
            nlcd_lc = read_nlcd_reclassified(case_name, paths["data"])

            modis_af, pxid_modis = rasterize_firms_detections(
                firms_modis, ig_date, lon, lat, nlcd_lc, cfg
            )
            viirs_af, pxid_viirs = rasterize_firms_detections(
                firms_viirs, ig_date, lon, lat, nlcd_lc, cfg
            )

            af_buffer = binary_dilation(
                (modis_af.astype(bool) | viirs_af.astype(bool)),
                structure=np.ones((cfg.af_buffer_size, cfg.af_buffer_size), dtype=bool),
            )
            af_contours = mask_to_contours(af_buffer)
            roi_bboxes = connected_component_bboxes(af_buffer)

            target_idx = select_center_roi_bbox(roi_bboxes, mtbs_bndy.shape)
            roi, roi_af, target_bbox = define_roi(
                mtbs_bndy, roi_bboxes, target_idx, cfg.roi_padding_pixels
            )
            vr = define_vr(
                roi, roi_af, nlcd_lc, target_bbox, cfg.vr_step_pixels, cfg.vr_max_steps
            )

            save_fire_map(
                paths["step1_figures"] / f"{case_name}.png",
                case_name, lon, lat, mtbs_bndy, pxid_modis, pxid_viirs,
                af_contours, mask_to_bbox(roi), mask_to_bbox(vr), cfg
            )

            save_step1_dataset(
                paths["step1_dataset"] / f"{case_name}.nc",
                case_name, mtbs_bndy, nlcd_lc, roi, vr, lon, lat, cfg
            )

            # -----------------------------------------------------------------
            # STEP 2. Preparation of analysis imagery
            # -----------------------------------------------------------------
            crop_bbox = mask_to_bbox(vr)
            lat, lon, mtbs_bndy, nlcd_lc, roi, vr = crop_base_layers(
                lat, lon, mtbs_bndy, nlcd_lc, roi, vr, crop_bbox
            )

            burn_boundaries = mask_to_contours(mtbs_bndy > 0)

            l8_si, l8_time = prepare_landsat_spectral_indices(
                case_name, pre_date, post_date, paths, crop_bbox,
                lon, lat, nlcd_lc, burn_boundaries, cfg
            )

            save_step2_dataset(
                paths["step2_dataset"] / f"{case_name}.nc",
                case_name, l8_si, l8_time, mtbs_bndy, nlcd_lc, roi, vr, lon, lat, cfg
            )

            # -----------------------------------------------------------------
            # STEP 3. Phenology detrending
            # -----------------------------------------------------------------
            l8_dsi, l8_pdsi, l8_time_pre, coef = run_phenology_detrending(
                case_name, ig_date, post_date, l8_si, l8_time, mtbs_bndy,
                nlcd_lc, vr, lon, lat, burn_boundaries, paths, cfg
            )

            save_step3_dataset(
                paths["step3_dataset"] / f"{case_name}.nc",
                case_name, l8_dsi, l8_pdsi, l8_time_pre, coef,
                mtbs_bndy, nlcd_lc, roi, vr, lon, lat, cfg
            )

            print(f"  -> Finished successfully: {case_name}")

        except Exception as exc:
            print(f'Case "{case_name}" failed.')
            print(str(exc))
            print(traceback.format_exc())


if __name__ == "__main__":
    main()
