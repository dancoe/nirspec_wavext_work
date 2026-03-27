#!/usr/bin/env python
"""
Determine Parlanti et al. 2025 Flux Correction Model Coefficients

This script solves for the throughput functions k(λ), a(λ), b(λ) using the Parlanti model:
    S(λ) ≈ k(λ)·f(λ) + a(λ)·f(λ/2) + b(λ)·f(λ/3)

Where:
    S(λ) = measured spectrum (M-grating with overlap contamination)
    f(λ) = intrinsic baseline spectrum (PRISM)
    k(λ) = 1st order throughput function
    a(λ) = 2nd order throughput function
    b(λ) = 3rd order throughput function

Usage:
    python determine_parlanti_coefficients.py
"""

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from stdatamodels.jwst import datamodels
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')


class ParlanthiCoefficientFitter:
    """Fit Parlanti model coefficients from spectral data."""
    
    def __init__(self, data_dir='/Users/dcoe/NIRSpec/wavext/data/PID1492',
                 output_dir='/Users/dcoe/NIRSpec/wavext/nirspec_wavext_work/plots'):
        self.data_dir = data_dir
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        
    def get_level3_spectrum(self, file_path):
        """Extract wavelength and flux from Level 3 (_x1d.fits) file."""
        try:
            model = datamodels.open(file_path)
            all_wav, all_flux = [], []
            specs = getattr(model, 'spec', [])
            if not specs and hasattr(model, 'spec_table'):
                specs = [model]
                
            for s in specs:
                wav = s.spec_table['wavelength']
                flx = s.spec_table['flux']
                # Filter valid values
                msk = (wav > 0) & (flx > 1e-12) & (~np.isnan(flx)) & (flx < 1e6)
                all_wav.append(wav[msk])
                all_flux.append(flx[msk])
            
            if len(all_wav) > 0:
                wav = np.concatenate(all_wav)
                flux = np.concatenate(all_flux)
                # Sort by wavelength
                idx = np.argsort(wav)
                return wav[idx], flux[idx]
        except Exception as e:
            print(f"Error reading {file_path}: {e}")
        return np.array([]), np.array([])
    
    def load_spectra(self):
        """Load PRISM and M-grating spectra."""
        print("Loading spectra...")
        
        # PRISM baseline
        prism_files = glob.glob(f'{self.data_dir}/*clear-prism*x1d.fits')
        if not prism_files:
            prism_files = glob.glob(f'{self.data_dir}/*0310e*x1d.fits')
        
        if prism_files:
            self.prism_wav, self.prism_flux = self.get_level3_spectrum(prism_files[0])
            print(f"  PRISM: {len(self.prism_wav)} points, λ=[{self.prism_wav[0]:.3f}, {self.prism_wav[-1]:.3f}] µm")
        else:
            raise FileNotFoundError("No PRISM data found")
        
        # M-grating spectra
        self.gratings = {}
        
        # Try multiple patterns for each grating
        grating_patterns = {
            'G140M': ['*g140m*extract_1d.fits', '*f100lp-g140m*x1d.fits', '*g140m*x1d.fits'],
            'G235M': ['*g235m*extract_1d.fits', '*f170lp-g235m*x1d.fits', '*g235m*x1d.fits'],
            'G395M': ['*g395m*extract_1d.fits', '*f290lp-g395m*x1d.fits', '*g395m*x1d.fits'],
        }
        
        for grating, patterns in grating_patterns.items():
            for pattern in patterns:
                files = glob.glob(f'{self.data_dir}/{pattern}')
                if files:
                    # Use the first file found
                    wav, flux = self.get_level3_spectrum(files[0])
                    if len(wav) > 0:
                        self.gratings[grating] = {'wav': wav, 'flux': flux, 'file': files[0]}
                        print(f"  {grating}: {len(wav)} points, λ=[{wav[0]:.3f}, {wav[-1]:.3f}] µm from {os.path.basename(files[0])}")
                        break
        
        if not self.gratings:
            raise ValueError("No M-grating data found")
    
    def create_interpolator(self, wav, flux, kind='linear', fill_value='extrapolate'):
        """Create interpolation function for spectrum."""
        return interp1d(wav, flux, kind=kind, bounds_error=False, fill_value=fill_value)
    
    def parlanti_model(self, wav, f_interp, k_coeffs, a_coeffs, b_coeffs):
        """
        Evaluate Parlanti model at given wavelengths.
        
        S(λ) = k(λ)·f(λ) + a(λ)·f(λ/2) + b(λ)·f(λ/3)
        
        where coefficients are provided as polynomials or interpolators.
        """
        # Get baseline at first order
        f1 = f_interp(wav)
        
        # Get baseline at second and third order wavelengths (shifted)
        f2 = f_interp(wav / 2.0)  # 2nd order at λ/2
        f3 = f_interp(wav / 3.0)  # 3rd order at λ/3
        
        # Evaluate coefficients (assumed to be polynomial)
        k_wav = np.polyval(k_coeffs, wav)
        a_wav = np.polyval(a_coeffs, wav)
        b_wav = np.polyval(b_coeffs, wav)
        
        return k_wav * f1 + a_wav * f2 + b_wav * f3
    
    def fit_coefficients_window(self, grating, wav_range=None, poly_order=3, plot=True):
        """
        Fit coefficients for a specific wavelength window using polynomial basis.
        
        Args:
            grating: Grating name (G140M, G235M, G395M)
            wav_range: Tuple (wav_min, wav_max) or None for full range
            poly_order: Polynomial order for fitting coefficients
            plot: Whether to generate diagnostic plots
        """
        if grating not in self.gratings:
            print(f"  {grating} not found")
            return None
        
        # Get M-grating spectrum
        s_wav = self.gratings[grating]['wav']
        s_flux = self.gratings[grating]['flux']
        
        # Create PRISM interpolator
        f_interp = self.create_interpolator(self.prism_wav, self.prism_flux)
        
        # Define wavelength window
        if wav_range is None:
            wav_range = (s_wav[0], s_wav[-1])
        
        # Select data in window
        mask = (s_wav >= wav_range[0]) & (s_wav <= wav_range[1])
        wav_window = s_wav[mask]
        flux_window = s_flux[mask]
        
        print(f"\n  Fitting {grating} in range λ=[{wav_range[0]:.3f}, {wav_range[1]:.3f}] µm")
        print(f"    Data points: {np.sum(mask)}")
        
        # Build design matrix for polynomial coefficients
        # We solve: S(λ) = k(λ)·f(λ) + a(λ)·f(λ/2) + b(λ)·f(λ/3)
        # Where k(λ), a(λ), b(λ) are polynomials of order poly_order
        
        n_coeffs = poly_order + 1
        n_total = 3 * n_coeffs  # k, a, b each have n_coeffs
        
        # Design matrix A where each row is [f(λ)·P(λ), f(λ/2)·P(λ), f(λ/3)·P(λ)]
        A = np.zeros((len(wav_window), n_total))
        
        f1 = f_interp(wav_window)
        f2 = f_interp(wav_window / 2.0)
        f3 = f_interp(wav_window / 3.0)
        
        # Polynomial basis for wavelength dependence
        for i in range(n_coeffs):
            power = poly_order - i
            poly_basis = wav_window ** power
            
            # Columns for k coefficients
            A[:, i] = f1 * poly_basis
            # Columns for a coefficients
            A[:, n_coeffs + i] = f2 * poly_basis
            # Columns for b coefficients
            A[:, 2*n_coeffs + i] = f3 * poly_basis
        
        # Solve least squares: A·x = S(λ)
        try:
            x, residuals, rank, s_values = np.linalg.lstsq(A, flux_window, rcond=None)
            
            # Extract coefficient polynomials
            k_coeffs = x[:n_coeffs]
            a_coeffs = x[n_coeffs:2*n_coeffs]
            b_coeffs = x[2*n_coeffs:3*n_coeffs]
            
            # Calculate fit quality
            s_model = self.parlanti_model(wav_window, f_interp, k_coeffs, a_coeffs, b_coeffs)
            rms_error = np.sqrt(np.mean((s_model - flux_window) ** 2))
            mean_flux = np.mean(flux_window)
            rms_error_pct = 100 * rms_error / mean_flux if mean_flux > 0 else 0
            
            print(f"    RMS error: {rms_error:.3e} ({rms_error_pct:.2f}%)")
            print(f"    Matrix rank: {rank}, singular values: {s_values[:3]}")
            
            result = {
                'grating': grating,
                'wav_range': wav_range,
                'k_coeffs': k_coeffs,
                'a_coeffs': a_coeffs,
                'b_coeffs': b_coeffs,
                'x_all': x,
                'residuals': residuals,
                'rank': rank,
                'rms_error': rms_error,
                'rms_error_pct': rms_error_pct,
                'wav_window': wav_window,
                'flux_window': flux_window,
                's_model': s_model,
            }
            
            if plot:
                self._plot_fit(result, f_interp)
            
            return result
            
        except np.linalg.LinAlgError as e:
            print(f"    ERROR: Linear algebra failed: {e}")
            return None
    
    def _plot_fit(self, result, f_interp):
        """Generate diagnostic plot of fit."""
        grating = result['grating']
        wav_window = result['wav_window']
        flux_window = result['flux_window']
        s_model = result['s_model']
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        fig.suptitle(f'{grating} Parlanti Model Fit\nλ=[{result["wav_range"][0]:.2f}, {result["wav_range"][1]:.2f}] µm', 
                     fontsize=14, fontweight='bold')
        
        # Plot 1: Data vs Model
        ax = axes[0, 0]
        ax.plot(wav_window, flux_window, 'k-', linewidth=1.5, label='Measured S(λ)', alpha=0.7)
        ax.plot(wav_window, s_model, 'r--', linewidth=1, label='Model fit', alpha=0.7)
        ax.set_xlabel('Wavelength (µm)')
        ax.set_ylabel('Flux (Jy)')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(alpha=0.3)
        ax.set_title('Spectrum vs Model')
        
        # Plot 2: Residuals
        ax = axes[0, 1]
        residuals = flux_window - s_model
        ax.plot(wav_window, residuals, 'b-', linewidth=0.5)
        ax.axhline(0, color='r', linestyle='--', alpha=0.5)
        ax.set_xlabel('Wavelength (µm)')
        ax.set_ylabel('Residuals (Jy)')
        ax.grid(alpha=0.3)
        ax.set_title(f'Residuals (RMS: {result["rms_error"]:.2e}, {result["rms_error_pct"]:.1f}%)')
        
        # Plot 3: Coefficient functions
        ax = axes[1, 0]
        k = result['k_coeffs']
        a = result['a_coeffs']
        b = result['b_coeffs']
        
        # Evaluate coefficients over range
        wav_eval = np.linspace(result['wav_range'][0], result['wav_range'][1], 100)
        k_wav = np.polyval(k, wav_eval)
        a_wav = np.polyval(a, wav_eval)
        b_wav = np.polyval(b, wav_eval)
        
        ax.plot(wav_eval, k_wav, 'b-', label='k(λ) [1st order]', linewidth=2)
        ax.plot(wav_eval, a_wav, 'g-', label='a(λ) [2nd order]', linewidth=2)
        ax.plot(wav_eval, b_wav, 'r-', label='b(λ) [3rd order]', linewidth=2)
        ax.set_xlabel('Wavelength (µm)')
        ax.set_ylabel('Throughput')
        ax.axhline(0, color='k', linestyle='-', alpha=0.2, linewidth=0.5)
        ax.axhline(1, color='k', linestyle='--', alpha=0.2, linewidth=0.5)
        ax.legend()
        ax.grid(alpha=0.3)
        ax.set_title('Coefficient Functions')
        
        # Plot 4: Component contributions
        ax = axes[1, 1]
        f1 = f_interp(wav_window)
        f2 = f_interp(wav_window / 2.0)
        f3 = f_interp(wav_window / 3.0)
        
        k_wav = np.polyval(k, wav_window)
        a_wav = np.polyval(a, wav_window)
        b_wav = np.polyval(b, wav_window)
        
        contrib_1 = k_wav * f1
        contrib_2 = a_wav * f2
        contrib_3 = b_wav * f3
        
        ax.plot(wav_window, contrib_1, 'b-', label='1st order: k(λ)·f(λ)', linewidth=1)
        ax.plot(wav_window, contrib_2, 'g-', label='2nd order: a(λ)·f(λ/2)', linewidth=1)
        ax.plot(wav_window, contrib_3, 'r-', label='3rd order: b(λ)·f(λ/3)', linewidth=1)
        ax.set_xlabel('Wavelength (µm)')
        ax.set_ylabel('Flux (Jy)')
        ax.set_yscale('log')
        ax.legend()
        ax.grid(alpha=0.3)
        ax.set_title('Component Contributions')
        
        plt.tight_layout()
        output_file = os.path.join(self.output_dir, f'parlanti_fit_{grating}.png')
        plt.savefig(output_file, dpi=150)
        print(f"    Saved plot: {output_file}")
        plt.close()
    
    def save_coefficients(self, results):
        """Save fitted coefficients to a file."""
        output_file = os.path.join(self.output_dir, 'parlanti_coefficients.txt')
        
        with open(output_file, 'w') as f:
            f.write("# Parlanti et al. 2025 Flux Correction Model Coefficients\n")
            f.write("# Model: S(λ) = k(λ)·f(λ) + a(λ)·f(λ/2) + b(λ)·f(λ/3)\n")
            f.write("# where f(λ) is PRISM baseline, k/a/b are polynomial functions\n\n")
            
            for result in results:
                if result is None:
                    continue
                    
                f.write(f"\n{'='*60}\n")
                f.write(f"Grating: {result['grating']}\n")
                f.write(f"Wavelength Range: {result['wav_range'][0]:.4f} - {result['wav_range'][1]:.4f} µm\n")
                f.write(f"RMS Error: {result['rms_error']:.3e} ({result['rms_error_pct']:.2f}%)\n")
                f.write(f"{'='*60}\n\n")
                
                f.write("Coefficients (highest to lowest power):\n\n")
                
                f.write("k(λ) [1st order throughput]:\n")
                for i, coeff in enumerate(result['k_coeffs']):
                    power = len(result['k_coeffs']) - 1 - i
                    f.write(f"  λ^{power}: {coeff:+.6e}\n")
                
                f.write("\na(λ) [2nd order throughput]:\n")
                for i, coeff in enumerate(result['a_coeffs']):
                    power = len(result['a_coeffs']) - 1 - i
                    f.write(f"  λ^{power}: {coeff:+.6e}\n")
                
                f.write("\nb(λ) [3rd order throughput]:\n")
                for i, coeff in enumerate(result['b_coeffs']):
                    power = len(result['b_coeffs']) - 1 - i
                    f.write(f"  λ^{power}: {coeff:+.6e}\n")
                
                f.write("\n")
        
        print(f"\nCoefficients saved to: {output_file}")
    
    def run(self):
        """Run complete analysis."""
        print("\n" + "="*70)
        print("PARLANTI MODEL COEFFICIENT DETERMINATION")
        print("="*70)
        
        # Load spectra
        self.load_spectra()
        
        # Fit coefficients for each M-grating
        results = []
        
        # Define wavelength windows based on expected overlap regions
        windows = {
            'G140M': (1.8, 3.3),      # Overlap region for G140M
            'G235M': (3.0, 5.3),      # Overlap region for G235M
            'G395M': (4.2, 5.6),      # Overlap region for G395M
        }
        
        for grating in ['G140M', 'G235M', 'G395M']:
            if grating in self.gratings:
                wav_range = windows.get(grating)
                result = self.fit_coefficients_window(grating, wav_range=wav_range, poly_order=3)
                results.append(result)
        
        # Save results
        self.save_coefficients(results)
        
        print("\n" + "="*70)
        print("Analysis complete!")
        print("="*70 + "\n")


if __name__ == '__main__':
    fitter = ParlanthiCoefficientFitter()
    fitter.run()
