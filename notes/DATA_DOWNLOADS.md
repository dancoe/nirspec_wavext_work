# DATA_DOWNLOADS.md - NIRSpec Wavelength Extension Data Acquisition

This document outlines the data acquisition strategy for the NIRSpec Wavelength Extension project, the scripts used for downloading, and how to verify that all data is correctly present on disk.

## Download Strategy

The project utilizes `astroquery.mast` to retrieve `_rate.fits` files directly from the MAST archive. We prioritize large (>50MB) files to capture science exposures while skipping target acquisition (TA) and confirmation images.

### Target PIDs

| Category | PIDs | Mode | Targets |
| :--- | :--- | :--- | :--- |
| **v5 Calibration** | 6644, 6645 | SLIT, IFU | NGC2506-G31, P330E-C3 |
| **Science Validation** | 2654, 2186 | IFU | SDSSJ0749, SDSSJ0841, UGC-5101 |
| **Baseline Cal** | 1536-1538 | IFU | J1743045, G191-B2B, P330E |
| **Extended FS** | 1492 | SLIT | IRAS-05248-7007 |

## Acquisition Scripts

The following scripts in `analysis/acquisition/` handle data retrieval:

1.  **`download_v5_data.py`** (Current): Downloads science validation targets and the v5 calibration sources (6644/6645).
2.  **`download_ifu_rates.py`**: Robust downloader for standard star programs (1536-1538, 6645). Features 3-attempt retry and file size filtering.
3.  **`download_fs_data.py`**: Legacy script for downloading Fixed Slit data for PID 1492.
4.  **`query_ifu_mast.py`**: Discovery script for identifying IFU observations.

## Verifying Progress

To check if downloads are complete, use the **`audit_downloads.py`** script. It queries the MAST archive for the expected number of science rate files and compares them against unique filenames on your local disk.

### Running the Audit

```bash
python3 analysis/acquisition/audit_downloads.py
```

### Manual Check

You can also perform a quick manual check using these terminal commands:

```bash
# Count rate files in science target folders
ls -l data/PID2654*/*rate.fits | wc -l
ls -l data/PID2186*/*rate.fits | wc -l

# Check for 0-byte or partial downloads
find data/ -size 0 -name "*rate.fits"
```

## Directory Structure

Data is organized by PID and target name under the `data/` directory.

-   **`data/PID####_TARGETNAME/`**: Standard storage for most programs (e.g., `data/PID2654_SDSSJ0749`).
-   **`data/IFU/`**: Historically used for standard star IFU calibrations (e.g., `data/IFU/PID1537_G191-B2B`).

For more details on specific observation IDs and exposure times, see [DATA.md](DATA.md). For a full index of analysis scripts, see [INDEX.md](INDEX.md).
