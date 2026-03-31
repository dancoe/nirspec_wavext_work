# Running the JWST pipeline for NIRSpec data

# Environment
micromamba activate jwst_1.20.2

# PYTHONPATH
```export PYTHONPATH=/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext```

# CRDS Reference Files

```
CRDS_SERVER_URL = https://jwst-crds.stsci.edu
CRDS_PATH = ~/crds_cache
```

# Online Notebooks
Best practices for running the JWST pipeline for NIRSpec:

analysis/reduction/
* JWPipeNB-NIRSpec-IFU.ipynb
* JWPipeNB-NIRSpec-FS.ipynb

via https://github.com/spacetelescope/jwst-pipeline-notebooks/tree/main/notebooks/NIRSPEC
* MOS: https://github.com/spacetelescope/jwst-pipeline-notebooks/blob/main/notebooks/NIRSPEC/MOS/JWPipeNB-NIRSpec-MOS.ipynb
* FS: https://github.com/spacetelescope/jwst-pipeline-notebooks/blob/main/notebooks/NIRSPEC/FSlit/JWPipeNB-NIRSpec-FS.ipynb
* IFU: https://github.com/spacetelescope/jwst-pipeline-notebooks/blob/main/notebooks/NIRSPEC/IFU/JWPipeNB-NIRSpec-IFU.ipynb
* BOTS: https://github.com/spacetelescope/jwst-pipeline-notebooks/tree/main/notebooks/NIRSPEC/BOTS
