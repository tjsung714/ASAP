# ASAP burn severity analysis pipeline

This repository provides two implementations of the ASAP preprocessing workflow for burn severity analysis.

- **MATLAB version** saves step outputs as `.mat`
- **Python version** saves step outputs as `.nc` (NetCDF via `xarray`)

The workflow integrates three steps into a single pipeline.

1. Definition of the analysis area
2. Preparation of analysis imagery
3. Phenology detrending

## Recommended repository structure

Place the scripts in the project root so that they can resolve paths automatically.

```text
project_root/
├── asap_pipeline_reviewed.m
├── asap_pipeline_python.py
├── requirements_python.txt
├── requirements_matlab.txt
├── 00_Data/
│   ├── metadata.xlsx
│   ├── FIRMS/
│   │   ├── fire_M-C61.csv
│   │   └── fire_V-C2.csv
│   ├── Landsat/
│   │   └── <case_name>*.tif
│   ├── NLCD/
│   │   └── <case_name>_*.tif
│   └── MTBS_WFIGS/
│       └── <case_name>
├── 01_ASAP_step1/
├── 02_ASAP_step2/
└── 03_ASAP_step3/
```

The output folders are created automatically if they do not already exist.

---

## Python version

### Environment requirements

- Python **3.10 or later**
- Packages listed in `requirements_python.txt`

### Install dependencies

Create and activate a virtual environment.

**Linux / macOS**

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements_python.txt
```

**Windows (PowerShell)**

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements_python.txt
```

### Run the pipeline

From the project root:

```bash
python asap_pipeline_python.py
```

### Python outputs

The Python pipeline saves one NetCDF file per wildfire case for each step.

```text
01_ASAP_step1/Dataset/<case_name>.nc
02_ASAP_step2/Dataset/<case_name>.nc
03_ASAP_step3/Dataset/<case_name>.nc
```

Figures are saved as PNG files under each step directory.

### Read the outputs in Python

Example for step 2:

```python
import xarray as xr

ds = xr.open_dataset("02_ASAP_step2/Dataset/<case_name>.nc")
print(ds)
ndvi = ds["NDVI"]
```

Example for step 3:

```python
import xarray as xr

ds = xr.open_dataset("03_ASAP_step3/Dataset/<case_name>.nc")
print(ds["dSI"])
print(ds["pdSI"])
print(ds["COEF"])
```

### Notes for Python users

- The script assumes that `metadata.xlsx` contains at least the following columns:
  - `Incid_Name`
  - `Ig_Date`
  - `Post_Date`
  - `Pre_Date(const)`
- The Python version uses `NetCDF` instead of `.mat` because labeled dimensions and time coordinates are easier to handle in Python.
- If `netCDF4` is installed, compressed NetCDF output is used automatically.

---

## MATLAB version

### Environment requirements

See `requirements_matlab.txt` for the required MATLAB products.

At minimum, you should have:

- MATLAB
- Image Processing Toolbox
- Mapping Toolbox

### Run the pipeline

1. Place `asap_pipeline_reviewed.m` in the project root.
2. Open MATLAB.
3. Change the current folder to the project root.
4. Run the script:

```matlab
asap_pipeline_reviewed
```

If the script is not placed in the project root, edit `path_root` at the top of the script.

### MATLAB outputs

The MATLAB pipeline saves one MAT file per wildfire case for each step.

```text
01_ASAP_step1/Dataset/<case_name>.mat
02_ASAP_step2/Dataset/<case_name>.mat
03_ASAP_step3/Dataset/<case_name>.mat
```

Figures are saved as PNG files under each step directory.

### Notes for MATLAB users

- The MATLAB version is best suited for users who want to keep the original MATLAB-based workflow.
- The output format is `.mat`, which is convenient inside MATLAB but less convenient than NetCDF for labeled multi-dimensional analysis in Python.
- If you need Python-friendly outputs, use the Python version.

---

## Main differences between the MATLAB and Python versions

| Item | MATLAB | Python |
|---|---|---|
| Main script | `asap_pipeline_reviewed.m` | `asap_pipeline_python.py` |
| Output dataset format | `.mat` | `.nc` |
| Figure format | `.png` | `.png` |
| Best for | Existing MATLAB workflow | Python analysis and downstream processing |
| Array handling | Numeric arrays in MAT files | Labeled arrays via `xarray` |

---

## Troubleshooting

### Python

- **`ModuleNotFoundError`**: install dependencies with `pip install -r requirements_python.txt`
- **Excel reading error**: make sure `openpyxl` is installed
- **NetCDF write issue**: install `netCDF4`
- **GeoTIFF read issue**: make sure `rasterio` installed successfully in your environment

### MATLAB

- **`Undefined function 'imdilate'` or `bwboundaries`**: Image Processing Toolbox is missing
- **`Undefined function 'geotiffread'` or `geotiffinfo`**: Mapping Toolbox is missing
- **Path-related errors**: check whether the script is placed in the project root or edit `path_root`

---

## Suggested citation / archive setup

A common release pattern is:

- Code on **GitHub**
- Archived release on **Zenodo**

For reproducibility, archive:

- the script version used for the paper
- this README
- dependency files
- a short description of the expected input folder structure
