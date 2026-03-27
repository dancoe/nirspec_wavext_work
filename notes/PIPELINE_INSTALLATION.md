# JWST Pipeline Installation Instructions

These notes cover how to install the latest version of the JWST calibration pipeline using `micromamba`.

## Latest Stable Release: 1.20.2 (DMS Build B12.1.1)
*Released: Oct 31, 2025*

### Standard Installation
The recommended way to install a specific version is into a fresh environment with Python 3.13.

```bash
# 1. Create a new environment
micromamba create -n jwst_1.20.2 python=3.13

# 2. Activate the environment
micromamba activate jwst_1.20.2

# 3. Install the specific pipeline version
pip install jwst==1.20.2
```

---

## Important "Gotchas"

### 1. Python Compatibility
*   **Supported:** 3.11, 3.12, 3.13
*   **Unsupported:** Python 3.14 (pipeline is not yet compatible).

### 2. Environment Location
*   Environments may be split between `/Users/dcoe/miniconda3/envs/` and `/Users/dcoe/micromamba/envs/`.
*   If you can't find an environment, check both locations or use `micromamba env list --json`.

### 3. Dependency Version Traps
*   **Versions 1.15.1 – 1.16.1:** These can pull an incompatible version of `gwcs`.
*   *Fix:* Run `pip install 'gwcs<0.22'` to downgrade.

### 4. CRDS Setup (Outside STScI Network)
To avoid excessive re-downloads and slow performance, always ensure your CRDS environment variables are set. Add these to your `.zshrc` or run them in your active session:

```bash
export CRDS_PATH=$HOME/crds_cache
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
```

### 5. MacOS Compatibility
*   Installation on **MacOS Mojave 10.14** will fail due to `opencv-python` build issues.
*   JWST requires a C compiler (e.g., Xcode Command Line Tools) for dependencies.

---

## Verifying Installation
Once installed, verify the version and dependencies:
```bash
pip show jwst
# To list all versions in the environment:
# micromamba list
```
