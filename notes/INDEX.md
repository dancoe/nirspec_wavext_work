# NIRSpec Wavelength Extension Project Index

This file provides a summary and links to the various project files, sorted by their last modification time (most recent first).

## 📝 Notes (`notes/`)
Key documentation, logs, and planning files.

| Last Modified | File | Summary |
| :--- | :--- | :--- |
| 2026-03-29 00:28:55 | [INDEX.md](INDEX.md) | This project index and summary file. |
| 2026-03-29 00:24:44 | [OBSERVATIONS.md](OBSERVATIONS.md) | Log of processed observations and results. |
| 2026-03-29 00:05:56 | [CALIBRATION.md](CALIBRATION.md) | Details of the calibration process for wavelength extension. |
| 2026-03-29 00:05:54 | [PARLANTI.md](PARLANTI.md) | Overview of the Parlanti et al. ghost contamination model. |
| 2026-03-29 00:04:02 | [CHANGE_LOG.md](CHANGE_LOG.md) | Chronological log of major updates and fixes. |
| 2026-03-29 00:04:02 | [LATEST_WORK.md](LATEST_WORK.md) | Current status and immediate next steps. |
| 2026-03-28 23:58:54 | [INSTRUCTIONS.md](INSTRUCTIONS.md) | Step-by-step guide to running the customized pipeline. |
| 2026-03-28 23:53:52 | [PROMPT_LOG.md](PROMPT_LOG.md) | Record of AI assistant interaction history. |
| 2026-03-28 22:34:40 | [PARLANTI_IFU.md](PARLANTI_IFU.md) | Deep-dive into Parlanti et al. model for IFU data. |
| 2026-03-28 22:14:09 | [PARLANTI_PLOTS.md](PARLANTI_PLOTS.md) | Plotting specifications (colors, layouts) for Parlanti analysis. |
| 2026-03-27 14:41:53 | [FLATS.md](FLATS.md) | Documentation for flat-field calibration steps. |
| 2026-03-27 14:41:16 | [FLATS_FILES.md](FLATS_FILES.md) | Inventory of flat-field reference files. |
| 2026-03-27 14:33:49 | [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) | Roadmap for wavelength extension feature. |
| 2026-03-27 14:32:19 | [DATA.md](DATA.md) | Information about datasets used in calibration. |
| 2026-03-27 10:36:02 | [FLATS_PLOTS.md](FLATS_PLOTS.md) | Index and descriptions of flat-field diagnostic plots. |
| 2026-03-27 00:02:23 | [NOTES.md](NOTES.md) | General brainstorming and project context. |
| 2026-03-27 00:01:16 | [FORK.md](FORK.md) | Notes on forking the official JWST pipeline. |
| 2026-03-26 23:11:21 | [PIPELINE_RUNNING.md](PIPELINE_RUNNING.md) | Commands for local execution of pipeline steps. |
| 2026-03-26 23:00:17 | [WAVELENGTH_RANGES.md](WAVELENGTH_RANGES.md) | Definitions of extended wavelength ranges by grating/filter. |
| 2026-03-26 22:55:27 | [WAVEXT.md](WAVEXT.md) | High-level overview of the 'wavext' feature. |
| 2026-03-26 22:35:04 | [PIPELINE_UPDATES.md](PIPELINE_UPDATES.md) | Notes on required changes for upstream integration. |
| 2026-03-26 20:44:24 | [CONTRIBUTING.md](CONTRIBUTING.md) | Guidelines for contributing to the project. |
| 2026-03-26 20:30:45 | [PIPELINE_INSTALLATION.md](PIPELINE_INSTALLATION.md) | Environment setup (micromamba/conda) instructions. |

## 📊 Analysis (`analysis/`)
Scripts and notebooks for calibration and visualization.

### Calibration (`analysis/calibration/`)
| Last Modified | File | Purpose |
| :--- | :--- | :--- |
| 2026-03-28 22:35:54 | [plot_parlanti_calib.py](../analysis/calibration/plot_parlanti_calib.py) | Master script for Parlanti calibration visualization. |
| 2026-03-28 08:19:00 | [solve_parlanti_calib.py](../analysis/calibration/solve_parlanti_calib.py) | Main solver for ghost contamination coefficients. |
| 2026-03-28 07:48:39 | [download_calspec.py](../analysis/calibration/download_calspec.py) | Fetches CalSpec standards for throughput fitting. |

### Plotting (`analysis/plotting/`)
| Last Modified | File | Purpose |
| :--- | :--- | :--- |
| 2026-03-28 23:53:52 | [plot_parlanti_components_multi.py](../analysis/plotting/plot_parlanti_components_multi.py) | Ghost component decomposition and multi-panel plotting. |
| 2026-03-28 09:50:31 | [generate_reference_star_plots.py](../analysis/plotting/generate_reference_star_extraction_plots.py) | Validation against star standard spectra. |
| 2026-03-27 16:50:39 | [plot_mos_extraction_regions.py](../analysis/plotting/plot_mos_extraction_regions.py) | Visualizes MOS extraction slits on detectors. |
| 2026-03-27 16:39:48 | [generate_final_plots.py](../analysis/plotting/generate_final_plots.py) | Final figures for project reports and documentation. |
| 2026-03-27 16:34:49 | [plot_final_coeffs.py](../analysis/plotting/plot_final_coeffs.py) | Visualizes the derived ghost coefficients. |
| 2026-03-27 15:17:15 | [plot_3source_summary.py](../analysis/plotting/plot_3source_summary.py) | Summary plots for 3-source calibration solvers. |
| 2026-03-27 14:32:57 | [batch_plot_new.py](../analysis/plotting/batch_plot_new.py) | Automation for large-scale plot generation. |
| 2026-03-27 14:03:10 | [inspect_spectra.py](../analysis/plotting/inspect_spectra.py) | Interactive spectral inspection tool. |
| 2026-03-27 11:59:53 | [reorganize_plots.py](../analysis/plotting/reorganize_plots.py) | Utility to sort plots into organized directory structures. |
| 2026-03-27 10:10:29 | [plot_flats_all_modes.py](../analysis/plotting/plot_flats_all_modes.py) | Comprehensive S-flat plane visualizations across modes. |

### Utils (`analysis/utils/`)
| Last Modified | File | Purpose |
| :--- | :--- | :--- |
| 2026-03-27 12:18:15 | [generate_sflat_tables_api.py](../analysis/utils/generate_sflat_tables_api.py) | API-driven s-flat table generation. |
| 2026-03-27 12:17:05 | [find_missing_sflats.py](../analysis/utils/find_missing_sflats.py) | Audits files for missing calibration data. |
| 2026-03-27 12:15:11 | [generate_sflat_table.py](../analysis/utils/generate_sflat_table.py) | Tabulates available S-flat planes for reference. |
| 2026-03-27 12:13:29 | [get_sflats.py](../analysis/utils/get_sflats.py) | Helper to retrieve specific S-flat files. |
| 2026-03-27 12:13:12 | [get_all_sflats.py](../analysis/utils/get_all_sflats.py) | Utility to gather all reference S-flats. |
| 2026-03-27 12:04:00 | [format_sflats.py](../analysis/utils/format_sflats.py) | Reformatting script for flat field data. |
| 2026-03-27 10:52:10 | [audit_all_flats.py](../analysis/utils/audit_all_flats.py) | Global audit of flat-field files. |
| 2026-03-27 09:33:30 | [inspect_waverange.py](../analysis/utils/inspect_waverange.py) | Diagnostic tool for checking ASDF contents. |
| 2026-03-27 00:13:20 | [find_centroid.py](../analysis/utils/find_centroid.py) | Centroid finding for spectral trace alignment. |
| 2026-03-27 00:11:29 | [view_ref_file.py](../analysis/utils/view_ref_file.py) | Quick viewer for JWST reference files. |
| 2026-03-26 23:17:54 | [create_waverange_extended.py](../analysis/utils/create_waverange_extended.py) | Generates the custom ASDF reference file. |

## 🖼️ Plots (`plots/`)
Diagnostic and final visualizations.

| Last Modified | Folder/File | Description |
| :--- | :--- | :--- |
| 2026-03-28 22:00:14 | [MOS_EXTRACTIONS.md](../plots/MOS_EXTRACTIONS.md) | Notes and results for MOS-specific extractions. |
| 2026-03-28 09:51:26 | **extraction/** | Spectral extraction plots (most recent batch). |
| 2026-03-28 07:58:07 | [ANALYSIS_RESULTS.md](../plots/ANALYSIS_RESULTS.md) | Summary of calibration results and residual analyses. |
| 2026-03-27 15:13:14 | **Parlanti/** | Ghost model calibration results and sub-model plots. |
| 2026-03-27 14:33:06 | jw01537...extract_1d.png | Latest individual spectral extraction visualization. |
| 2026-03-27 13:45:13 | parlanti_fit_G140M.png | Throughput fit for G140M grating. |
| 2026-03-27 13:05:05 | **flats/** | Flat-field plane and extension visualizations. |

## 📚 References (`references/`)
External data and documentation.

| Last Modified | Item | Description |
| :--- | :--- | :--- |
| 2026-03-27 14:13:15 | [NIRSpec_fluxcal_observations.xlsx](../references/NIRSpec_fluxcal_observations.xlsx) | Master observation database. |
| 2026-03-27 10:48:06 | **Parlanti/** | Scientific context and reference papers. |
| 2026-03-27 10:47:58 | **flats/** | S-flat format documentation and screenshots. |

## 📂 Reference Files (`reference_files/`)
| Last Modified | File | Purpose |
| :--- | :--- | :--- |
| 2026-03-26 23:18:03 | [wavelengthrange_extended.asdf](../reference_files/wavelengthrange_extended.asdf) | Custom reference file for the pipeline. |

## ⚙️ Scripts (`scripts/`)
Workflow automation and debugging tools.

| Last Modified | File | Purpose |
| :--- | :--- | :--- |
| 2026-03-27 23:33:08 | [debug_nrs1.py](../scripts/debug_nrs1.py) | Investigating NRS1-specific coordinate issues. |
| 2026-03-27 23:30:16 | [test_discovery.py](../scripts/test_discovery.py) | Validates pipeline local-file discovery logic. |
| 2026-03-27 23:28:28 | [generate_1492_extractions.py](../scripts/generate_1492_extractions.py) | Automated extraction generator for PID 1492. |
| 2026-03-27 22:13:21 | [generate_missing_extractions.py](../scripts/generate_missing_extractions.py) | Utility to process datasets that failed in initial runs. |
