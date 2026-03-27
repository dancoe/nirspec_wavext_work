# Contributing to JWST Calibration Pipeline (Distilled)

Source: https://github.com/spacetelescope/jwst/blob/main/CONTRIBUTING.md

## Reporting Issues
- **Bugs & Feature Requests:** [Open an issue on GitHub](https://github.com/spacetelescope/jwst/issues) or contact the [JWST Help Desk](https://jwsthelp.stsci.edu).

## Development Setup
1. **Fork and Clone:**
   - [Fork the repository](https://github.com/spacetelescope/jwst/fork).
   - Clone your fork locally:
     ```bash
     git clone https://github.com/YOUR_USERNAME/jwst
     cd jwst/
     ```
   - Add the upstream repository:
     ```bash
     git remote add upstream https://github.com/spacetelescope/jwst
     ```

2. **Environment:**
   - Create and activate a development environment:
     ```bash
     mamba create -n jwst_dev python=3.13
     mamba activate jwst_dev
     pip install -e .
     ```

3. **Code Quality Tools:**
   - Install `pre-commit` to automate formatting checks:
     ```bash
     pip install pre-commit
     pre-commit install
     ```
     (Manual run: `pre-commit run --all`)

## Contribution Workflow
1. **Branching:** Create a new branch for each feature or fix.
   ```bash
   git fetch upstream --tags
   git checkout upstream/main -b feature/your_feature_name
   ```

2. **Implementation:** Make your changes and commit them.
   ```bash
   git add <files>
   git commit -m "Description of changes"
   git push -u origin feature/your_feature_name
   ```

3. **Staying Current:** Rebase on `upstream/main` regularly.
   ```bash
   git fetch --all
   git rebase -i upstream/main
   git push -u origin -f feature/your_feature_name
   ```

4. **Pull Request:** Open a PR on GitHub. Complete the task checklist and ensure tests pass.

## Standards
- **Style:** PEP8 (enforced via `ruff`).
- **Documentation:** Numpy style docstrings.
- **Dependencies:** If changing a dependency (e.g., `stcal`), install it in editable mode too: `pip install -e ../stcal`.
