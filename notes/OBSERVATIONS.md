# NIRSpec Flux Calibration Observations

**Data Offload Note:** As of March 30, 2026, most Level 1 and Level 2 data (including `rate`, `cal`, and intermediate reduction stage files) have been offloaded to **central store** (`/Volumes/wit4/...`) to save local storage. 
*   **Local retention:** Level 3 products (`spec3_*`, `stage3*`) and final v7 Jy-calibrated products remain on the laptop for analysis.
*   **Central store:** Organised into `FS/` and `IFU/` subdirectories.
*   **IFU v7 note (2026-03-30):** For IFU mode, `photom_0016` has `relresponse=1.0` and `photmj=1.0` throughout — no photom extension needed (unlike FS). IFU v7 re-reductions are stage2+3 fresh runs using Parlanti reffiles for consistency. Rate files for 1536/1537/1538 re-downloaded from MAST; 2186 rate files were already local.

**Source:** `data/NIRSpec_fluxcal_observations.xlsx` (v1.1, 24 Dec 2025)  
**Generated:** 2026-03-29

This document catalogs all NIRSpec spectrophotometric calibration programs relevant to the
wavelength-extension project, based on the STScI fluxcal observation spreadsheet. It documents
modes, gratings, targets, spectral types, and prioritization for our FS and IFU analyses.

---

## Detailed FS Observation Matrix

| PID  | Target         | Gratings       | Slits            | Notes                                   |
|:-----|:---------------|:---------------|:-----------------|:----------------------------------------|
| 1492 | IRAS-05248     | G140M, G235M, G395M | S200A1, S400A1, S1600A1 | Hold-out validation AGN/Galaxy target. |
| 1536 | J1743045       | ALL M+H, PRISM | S1600A1          | Primary standard star A8III.            |
| 1537 | G191-B2B       | ALL M+H, PRISM | S1600A1          | Hot white dwarf standard WD-DA.8.       |
| 1538 | P330E          | ALL M+H, PRISM | S1600A1          | G-type G2V standard (break k/α degeneracy).|
| 1678 | J1757132, etc. | G235M, G395M   | S1600A1          | 4 A-stars, GO-CAL program.              |
| 6644 | NGC2506-G31    | G140M, G235M, G395M | S1600A1     | Essential cool G1V star for v7 calibration. |
| 2186 | UGC5101       | G235M, G395M        | IFU         | ULIRG science target (v9 validation cohort).|
| 2654 | SDSSJ0841       | G140M, G235M        | IFU         | Quasar science target (v8 validation cohort).|

---

| PID | Target(s) | FS L1 | FS L2 | FS L3 | FS v7* | IFU L1 | IFU L2 | IFU L3 | IFU v7* | Priority |
|:----|:----------|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---------|
| 1536 | J1743045 | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅¹ | ★★★ |
| 1537 | G191-B2B | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅¹ | ★★★ |
| 1538 | P330E | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅ | ✅¹ | ★★★ |
| 1492 | IRAS-05248 | ✅ | ✅ | ✅ | ✅ | — | — | — | — | ★★★ |
| 6644 | NGC2506-G31| ✅ | ✅ | ✅ | ✅ | — | — | — | — | ★★★★ |
| 2186 | UGC5101 | — | — | — | — | ✅ | ✅ | ✅ | ✅² | ★★★ |
| 6645 | P330E | — | — | — | — | 🔄 | 🔄 | 🔄 | — | ★★★ |
| 2654 | SDSSJ0749 | — | — | — | — | 🔄 | 🔄 | 🔄 | — | ★★★ |

*\*v7 indicates IFU stage2+3 re-reduction with Parlanti reffiles using pipeline v1.20.2.*  
*¹ Rate files re-downloaded from MAST (2026-03-30); stage2+3 queued/running.*  
*² Stage3 re-run from existing stage2_ext cal files (rate files already local).*

---

## Target Spectral Type Reference

| Target        | SpT     | CALSPEC? | Notes                                     |
|---------------|---------|----------|-------------------------------------------|
| BD+60-1753    | A1V     | Yes      | Monitoring star                           |
| C26202        | F8IV    |          |                                           |
| G191-B2B      | WD-DA.8 | Yes      | Hot white dwarf standard                  |
| GD153         | WD-DA1.2| Yes      | Hot white dwarf standard                  |
| GD71          | WD-DA1.5| Yes      | Hot white dwarf standard                  |
| HD163466      | A6V     |          | GO-CAL                                    |
| HD18511       | K0III   | No       | Cool star, not in CALSPEC                 |
| HD19467       | G3V     | No       | Not in CALSPEC (PSF star)                 |
| J1743045      | A8III   | Yes      | Gordon et al. standard                    |
| J1757132      | A3V     | Yes      | Gordon et al. standard                    |
| J1802271      | A2V     | Yes      | Gordon et al. standard                    |
| J1805292      | A4V     | Yes      | Gordon et al. standard                    |
| J1808347      | A3V     | Yes      | Gordon et al. standard                    |
| J1812095      | A3V     | Yes      | Gordon et al. standard                    |
| LDS749B       | WD-DBQ4 | Yes      | White dwarf standard                      |
| NGC2506-G31   | G1V     | Yes      | COOL G-type; key for k/α degeneracy break |
| P177D         | G0V     | Yes      | G-type standard                           |
| P330-E        | G2V     | Yes      | Gray/G-type standard; Parlanti primary    |
| PG-1057+719   | WD-DA1.2|          | Early checkout (unusable)                  |
| SF1615        | G       |          | G-type standard                           |
| SNAP-2        | G0-5    |          | Quasi-stellar standard analog             |
| WD1057+719    | WD-DA1.2|          |                                           |
| WDFS0122-30   | WD-DA   |          |                                           |
| WDFS0458-56   | WD-DA   |          |                                           |
| WDFS1304+10   | WD-DA   |          |                                           |
| WDFS1557+55   | WD-DA   |          |                                           |
| WDFS2317-29   | WD-DA   |          |                                           |

---

## Program Summary

### Cycle 0

#### PID 1118 (PI: Proffitt)
- **Target:** J1743045 (A8III), SNAP-2 (G0-5)
- **Modes:** FS (S200A1/S200A2), IFU
- **Gratings:** G140H/F100LP (FS); G140M/F100LP (IFU)
- **Notes:** Early Cycle 0 data; J1743045 FS + IFU G140M observations available

#### PID 1128 (PI: Lützgendorf)
- **Target:** J1808347 (A3V), SNAP-2 (G0-5), J1743045 (A8III), P177D (G0V), PG-1057+719 (WD)
- **Modes:** FS (S1600A1, S200A1, S200A2, S400A1), IFU, MOS
- **Gratings:** PRISM, G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP
- **Notes:** Very comprehensive Cycle 0 calibration; J1808347 observed in ALL gratings + ALL
  three FS modes + IFU + MOS. Observations 15–17 (P177D, PG-1057+719) marked as early
  checkout with temperature instability — **do not use**.
  J1808347 is an A3V CALSPEC star with full grating coverage, making it excellent for FS analysis.

---

### Cycle 1

#### PID 1487 (PI: Manjavacas)
- **Target:** J1743045 (A8III), SF1615 (G)
- **Modes:** FS (S200A1, S200A2, S400A1, S1600A1), IFU (SF1615)
- **Gratings:** G140M/F100LP, PRISM
- **Notes:** SF1615 IFU PRISM; J1743045 FS PRISM + G140M in multiple slits

#### PID 1493 (PI: Proffitt)
- **Target:** NGC2506-G31 (G1V)
- **Modes:** MOS
- **Gratings:** PRISM
- **Notes:** G1V cool star in MOS; limited for FS analysis

#### PID 1536 (PI: Gordon) ★★★ **PRIMARY**
- **Target:** J1743045 (A8III), J1802271 (A2V), J1757132 (A3V)
- **Modes:** FS (S1600A1) + IFU (J1743045 only)
- **Gratings:** PRISM, G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP  
- **Notes:** Gordon primary fluxcal program Cycle 1. **IFU already reduced (G140M+G235M+G395M)**.
  FS data downloaded (MAST L3 x1ds + per-exposure x1d + NRS2 rate files).
  **NRS2 rate files in hand for G140M/F100LP, G235M/F170LP, G395M/F290LP → ready to extract.**

#### PID 1537 (PI: Gordon) ★★★ **PRIMARY**
- **Target:** G191-B2B (WD-DA.8), GD153 (WD-DA1.2)
- **Modes:** FS (S1600A1) + IFU (G191-B2B)
- **Gratings:** PRISM, G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP
- **Notes:** Gordon primary fluxcal program Cycle 1. **IFU already reduced**. FS:
  MAST L3 x1ds exist + **3 NRS2 extract_1d files already extracted**:
  - G140M NRS2: dither 3 → 1.96–3.16 µm ✅
  - G235M NRS2: dither 3 → 3.28–5.31 µm ✅  
  - G395H NRS2: dither 2 → 5.44–5.62 µm (edge, low utility)

#### PID 1538 (PI: Gordon) ★★★ **PRIMARY**
- **Target:** P330-E (G2V), P177D (G0V), NGC2506-G31 (G1V)
- **Modes:** FS (S1600A1) + IFU (P330-E)
- **Gratings:** PRISM, G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP
- **Notes:** Gordon primary fluxcal Cycle 1. P330-E is a G2V star (different from hot WDs/A-stars)
  — **critically important for breaking k/α degeneracy**. **IFU already reduced**. FS:
  MAST L3 x1ds (Obs 60 = G140M+G235M+G395M) + **3 NRS2 extract_1d files already extracted**:
  - G140M NRS2 (Obs 160, dither 1) → 1.96–3.16 µm ✅
  - G235M NRS2 (Obs 160, dither 5) → 3.28–5.31 µm ✅
  - G395H NRS2 (Obs 160, dither 4) → 5.45–5.62 µm (edge)
  NGC2506-G31 Obs 32 PRISM only; P177D Obs 4 PRISM only.

#### PID 1539 (PI: Gordon)
- **Target:** BD+60-1753 (A1V), SF1615 (G)
- **Modes:** FS monitoring
- **Gratings:** G140H/F100LP
- **Notes:** Monitoring program; BD+60-1753 observed repeatedly as calibration check star.

#### PID 1678 (PI: Ashby) ★★
- **Target:** J1757132 (A3V), J1812095 (A3V), J1802271 (A2V), J1805292 (A4V)
- **Modes:** FS (S1600A1) only
- **Gratings:** G235M/F170LP, G395M/F290LP (Obs 5–8 each)
- **Notes:** GO-CAL Cycle 1 program targeting A-type stars. Four CALSPEC stars with
  G235M + G395M coverage — useful as additional calibrators for G235M NRS2 extension.
  All targets are A-stars (similar to our primary set), so limited for k/α degeneracy.
  **Not yet downloaded.**

---

### Cycle 2

#### PID 3399 (PI: Perrin) ★★
- **Target:** J1757132 (A3V), HD163466 (A6V)
- **Modes:** IFU
- **Gratings:** G140H/F100LP, G235H/F170LP, G395H/F290LP
- **Notes:** GO-CAL Cycle 2 IFU. HD163466 is A6V, CALSPEC? GO-CAL designates it as a
  calibration target. Two sources with full H/M coverage in IFU. Would add data points
  but still A-stars. **Not yet downloaded.**

#### PID 4460 / 4465 (PI: Manjavacas)
- **Target:** J1743045 (A8III)
- **Modes:** FS (S200A1, S200A2, S400A1)
- **Gratings:** G140M/F100LP
- **Notes:** Cycle 2 repeat of Cycle 1 obs 1487; FS G140M only.

#### PID 4496 (PI: Gordon)
- **Target:** LDS749B (WD-DBQ4)
- **Modes:** Various
- **Gratings:** PRISM
- **Notes:** White dwarf LDS749B in Cycle 2; PRISM only.

#### PID 4497 (PI: Gordon)
- **Target:** J1812095 (A3V)
- **Modes:** FS + IFU
- **Gratings:** PRISM
- **Notes:** Monitoring/repeat.

#### PID 4498 (PI: Gordon) ★★
- **Target:** P330-E (G2V)
- **Modes:** FS (S1600A1, ALLSLITS, FULL) + IFU
- **Gratings:** G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP
- **Notes:** P330-E Cycle 2 repeat of PIDs 1538 FS+IFU Gordon program. Provides
  temporal baseline check. **Not yet downloaded but high scientific value as G2V.**

#### PID 4499 (PI: Gordon)
- **Target:** BD+60-1753 (A1V)
- **Modes:** FS monitoring
- **Gratings:** G140H/F100LP
- **Notes:** Monitoring.

---

### Cycle 3

#### PID 6605 (PI: Gordon)
- **Target:** J1802271 (A2V)
- **Modes:** FS
- **Gratings:** Various (from Gordon Cycle 3 program)
- **Notes:** Repeat fluxcal.

#### PID 6606 (PI: Gordon) ★★★
- **Target:** P330-E (G2V)
- **Modes:** FS (S1600A1, ALLSLITS) + IFU
- **Gratings:** G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP
- **Notes:** P330-E Cycle 3 fluxcal. Multiple observations per grating (Obs 31=deep,
  Obs 32=ALLSLITS). Same target as PIDs 1538, 4498 — provides 3 epochs of G2V data.
  **Not yet downloaded.**

#### PID 6607 (PI: Gordon)
- **Target:** BD+60-1753 (A1V)
- **Modes:** FS monitoring
- **Gratings:** G140H/F100LP, G235H/F170LP, G395H/F290LP
- **Notes:** Monitoring.

#### PID 6643 (PI: Proffitt) ★
- **Target:** GD71 (WD-DA1.5), G191-B2B (WD-DA.8)
- **Modes:** MOS (Q4 Field Point)
- **Gratings:** G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP
- **Notes:** Cycle 3 MOS fluxcal with two WD standards. MOS mode is different from FS/IFU;
  useful to cross-check but requires different analysis path.

#### PID 6644 (PI: Proffitt) ★★★★ **HIGH PRIORITY — COOL STARS**
- **Target:** NGC2506-G31 (G1V), HD18511 (K0III)
- **Modes:** FS (S1600A1)
- **Gratings:** G140M/F070LP, G140M/F100LP, G235M/F170LP, G395M/F290LP
- **Notes:** **This is the most important program for breaking the k/α degeneracy!**
  - NGC2506-G31 (G1V) IS in CALSPEC — a cool G-type star with very different f(λ/2)/f(λ)
    ratio compared to A-stars/WDs. Adding it to the 3-source solver will allow separation
    of k(λ) from α̃(λ). 
  - HD18511 (K0III) is NOT in CALSPEC — limited utility for Parlanti solver but useful
    as a spectral shape diagnostic.
  - NGC2506-G31 is also in PID 1538 Obs 32 (PRISM only) and PID 1538 Obs 32.
  **Priority: download and reduce NGC2506-G31 FS data. Not yet downloaded.**

#### PID 6645 (PI: Zeidler) ★★★
- **Target:** P330-E (G2V)
- **Modes:** IFU
- **Gratings:** G140M/F070LP, G140M/F100LP, G140H/F070LP, G140H/F100LP,
  G235M/F170LP, G235H/F170LP, G395M/F290LP, G395H/F290LP
- **Notes:** P330-E Cycle 3 IFU, Parlanti's "4th source" program. Multiple pointing positions
  (Obs 1 = star at edge, Obs 3 = star centered). This data was used by Parlanti et al. as
  the 4th IFU calibrator. **Currently being reduced (PID 6645 in IFU pipeline queue).**

---

## Mode Coverage Matrix

### Fixed Slit (FS) Programs by Grating

| Grating/Filter | PIDs with FS S1600A1 (CALSPEC) |
|:---------------|:-------------------------------|
| G140M/F070LP   | 1536, 1538, 4498, 6606, 7615   |
| G140M/F100LP   | 1536, 1537, 1538, 4498, 6606, 7615 |
| G140H/F070LP   | 1536, 1537, 1538, 4498, 6606, 7615 |
| G140H/F100LP   | 1536, 1537, 1538, 4498, 6606   |
| G235M/F170LP   | 1536, 1537, 1538, 1678, 4498, 6606, 6644, 7565, 7615 |
| G235H/F170LP   | 1536, 1537, 1538, 4498, 6606   |
| G395M/F290LP   | 1536, 1537, 1538, 1678, 4498, 6606, 6644, 7565 |
| G395H/F290LP   | 1536, 1537, 1538, 4498, 6606   |

### IFU Programs by Grating

| Grating/Filter | PIDs with IFU (CALSPEC)         |
|:---------------|:--------------------------------|
| G140M/F070LP   | 1128, 6606, 6645                |
| G140M/F100LP   | 1118, 1536, 1537, 1538, 4498, 6606, 6645 |
| G140H/F100LP   | 1536, 1537, 1538, 3399, 4498, 6606 |
| G235M/F170LP   | 1536, 1537, 1538, 4498, 6606, 6645 |
| G235H/F170LP   | 1536, 1537, 1538, 3399, 4498, 6606, 6645 |
| G395M/F290LP   | 1536, 1537, 1538, 4498, 6606, 6645 |
| G395H/F290LP   | 1536, 1537, 1538, 3399, 4498, 6606, 6645 |

---

## Priority Analysis for This Work

### Immediate Priority (data in hand)

| Source | PID | Mode | Status |
|:-------|:----|:-----|:-------|
| G191-B2B | 1537 | FS | ✅ v7 Jy extracted |
| P330E | 1538 | FS | ✅ v7 Jy extracted |
| J1743045 | 1536 | FS | ✅ v7 Jy extracted |
| NGC2506-G31| 6644 | FS | ✅ v7 Jy extracted (G1V) |
| IRAS-05248 | 1492 | FS | ✅ v7 Jy extracted (Hold-out Val) |
| UGC5101 | 2186 | IFU (G235M, G395M) | 🔄 v7 reduction pending |

The 5 FS sources (G191-B2B, P330E, J1743045, NGC2506-G31, and 1492) are now fully reprocessed through the v7 pipeline with the extended photom reference files, providing the foundation for the v7 FS calibration.

### Second Priority (download needed)

| Target | PID | SpT | Key Value |
|:-------|:----|:----|:----------|
| NGC2506-G31 | 6644 | G1V | Break k/α degeneracy — CALSPEC cool star |
| NGC2506-G31 | 1538 Obs 32 | G1V | PRISM only — limited but available |
| P330-E | 6606 | G2V | Cycle 3 repeat, temporal check |
| J1808347 | 1128 | A3V | Comprehensive grating coverage |
| G191-B2B | 7565 | WD | Cycle 4 FS S200A1/A2/A1 — multiple slits |

### Third Priority (A-star multiplicity)

| Program | Targets | Value |
|:--------|:--------|:------|
| PID 1678 | J1757132, J1812095, J1802271, J1805292 | 4 A-stars, G235M+G395M only |
| PID 3399 | J1757132, HD163466 | 2 A-stars IFU, G140H+G235H+G395H |
| PID 4498 | P330-E | G2V, Cycle 2 FS+IFU all gratings |

---

## Scientific Rationale for Prioritization

### Why Cool Stars Matter

The Parlanti k/α degeneracy arises because all three primary calibrators (G191-B2B, 
P330-E, J1743045) are dominated by hot or warm continua where f(λ/2)/f(λ) has nearly
the same value for all three (power-law-like SEDs). This means the NNLS solver cannot 
distinguish whether the extended-range flux deficit comes from k < 1 (first-order 
throughput) or α > 0 (second-order contamination).

Adding NGC2506-G31 (G1V) breaks this because:
1. G1V stars have **molecular absorption at λ/2** (CO bands at ~2.3+ µm → showing 
   up at λ~4.6 µm in second order). This makes f(λ/2)/f(λ) strongly wavelength-dependent
   in a uniquely different way from hot stars.
2. The G1V continuum has a steeper Rayleigh-Jeans slope than hot stars, so the ratio 
   f(λ/2)/f(λ) is different in amount as well.
3. This provides genuinely independent equations to uniquely solve for k, α, β.

### AGN Targets (from Parlanti et al.)

Parlanti et al. used two AGN programs for validation:
- **PID 2654**: SDSSJ0749, SDSSJ0841 (for G140M validation)
- **PID 2186**: UGC-5101 and 2 other local ULIRGs (for G235M validation)

These targets are not flux standards, but their spectra in the nominal range (NRS1) serve
as self-calibration truth. If future AGN data becomes available (e.g., from our PID 1492
IRAS-05248), it can be used as a cross-check after the k/α/β calibration is established
from standard stars.

Additionally, bright AGN observed in consecutive gratings (G140M + G235M or G235M + G395M)
make excellent validation targets. We should keep an eye out for any such datasets.

---
