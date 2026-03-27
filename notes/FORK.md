# JWST Pipeline Fork Details

To contribute to the JWST Pipeline, follow these steps to manage your local development environment.

## 1. Fork on GitHub
You have successfully forked the official repository to:
- **Fork URL:** [https://github.com/dancoe/jwst_nirspec_wavext](https://github.com/dancoe/jwst_nirspec_wavext)
- **Account:** dancoe@gmail.com (GitHub user: `dancoe`)

## 2. Local Setup (Command Line)
The repository has been successfully cloned using HTTPS to:
`~/NIRSpec/wavext/jwst_nirspec_wavext`

### Setup Commands used:
```bash
cd ~/NIRSpec/wavext/
git clone https://github.com/dancoe/jwst_nirspec_wavext.git
cd jwst_nirspec_wavext/
git remote add upstream https://github.com/spacetelescope/jwst
git fetch upstream
git checkout -b feature/nirspec_wavelength_extension
```

## 3. Pushing Changes
- You will be pushing code to **your fork** on GitHub (`origin`).
- From GitHub, you will then open a **Pull Request** to merge your changes into the official `spacetelescope/jwst` repository (`upstream`).

```bash
git add <modified_files>
git commit -m "Brief description of updates"
git push origin feature/nirspec_wavelength_extension
```
