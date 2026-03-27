# JWST Pipeline Development Notes

This repository (`nirspec_wavext_work`) contains key documentation, plans, analysis, and reference files for the **NIRSpec Wavelength Extension ("wavext")** project.

## Directory Structure
- **[notes/](notes/):** Project documentation and logs (this directory).
- **[analysis/](analysis/):** Data processing and calibration scripts.
- **[reference_files/](reference_files/):** Custom ASDF reference files for the extended pipeline.

## Project Management (notes/)
- **[INSTRUCTIONS.md](notes/INSTRUCTIONS.md):** Persistent instructions for maintaining project logs and development standards.
- **[PROMPT_LOG.md](notes/PROMPT_LOG.md):** Numbered history of all user prompts and instructions.
- **[CHANGE_LOG.md](notes/CHANGE_LOG.md):** Detailed log of all code and configuration changes, newest first.

## Setup & Implementation (notes/)
- **[PIPELINE_INSTALLATION.md](notes/PIPELINE_INSTALLATION.md):** Instructions for setting up and managing JWST pipeline environments.
- **[FORK.md](notes/FORK.md):** Details about the `jwst_nirspec_wavext` fork and local repository configuration.
- **[CONTRIBUTING.md](notes/CONTRIBUTING.md):** Distilled version of the official STScI JWST contributing guide.

## Research & Strategy (notes/)
- **[WAVEXT.md](notes/WAVEXT.md):** General project plan, resource links, and research summaries for the extension.
- **[PIPELINE_UPDATES.md](notes/PIPELINE_UPDATES.md):** Technical plan and notes for removing hard-coded limitations in `assign_wcs` for NRS2.
- **[WAVELENGTH_RANGES.md](notes/WAVELENGTH_RANGES.md):** Reference for adopted wavelength ranges from `msaexp` and Parlanti et al.
- **[DATA.md](notes/DATA.md):** Summary of observations for calibration (e.g., PID 1492).
- **[CALIBRATION.md](notes/CALIBRATION.md):** Mathematical models and strategy for throughput solve.
