# JWST Development Instructions

We've forked the JWST pipeline:  
`~/NIRSpec/wavext/jwst_nirspec_wavext`.

And we're adding NIRSpec a wavelength extension feature:
`feature/nirspec_wavelength_extension`

Please follow these instructions:

- **Current Status:** Refer to [LATEST_WORK.md](LATEST_WORK.md) for our most recent progress and upcoming plans; for more details on the long-term project goals and technical implementation, see [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md).
- **Log all prompts:** Record every user request numerically with a timestamp in [PROMPT_LOG.md](PROMPT_LOG.md).
- **Log all changes:** Detail every modification to the codebase in [CHANGE_LOG.md](CHANGE_LOG.md) (newest on top).
- **Maintain documentation:** Keep these instructions and other notes up to date.
- **Commit changes:** After making updates, commit and push changes to the repositories (`nirspec_wavext_work`) and/or (`jwst_nirspec_wavext`) as warranted.
- **Don't operate the browser** unless I ask you to. The terminal should suffice.
## Pipeline Environment
All wavelength extension processing MUST be performed in the `jwst_1.20.2` environment. This version is confirmed compatible with the current codebase and supports the necessary `stdatamodels` and `asdf` schema resolutions.

```bash
# Activation
micromamba activate jwst_1.20.2

# Using custom patch
export PYTHONPATH=/Users/dcoe/NIRSpec/wavext/jwst_nirspec_wavext
```

## Processing Workflow
1.  **Preparation:** Download Level 3 (f(λ)) and Level 2 (obs) data for the target (e.g., PID 1492, 1537, 1538).
2.  **Assignment:** Run `assign_wcs` with the `wavelengthrange_extended.asdf` override.
3.  **Extraction:** Run `extract_2d` and `extract_1d` (to generate extended spectral cutouts).
4.  **Analysis:** Use the baseline f(λ) from nominal extractions (or PRISM) to solve for the throughput functions k(λ), α̃(λ), β̃(λ) using the multi-source Parlanti model.
