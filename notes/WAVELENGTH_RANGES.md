# Wavelength Ranges

## Adopted Ranges for WAVEXT

These ranges will be used to populate the custom `wavelengthrange` reference file.

| Filter_Grating | Nominal (µm) | Proposed (µm) | Notes |
| :--- | :--- | :--- | :--- |
| `F070LP_G140M/H` | 0.7 – 1.8 | 0.6 – 3.3 | Extended for NRS2 coverage |
| `F100LP_G140M/H` | 1.0 – 1.8 | 0.9 – 3.3 | Extended for NRS2 coverage |
| `F170LP_G235M/H` | 1.7 – 3.1 | 1.5 – 5.3 | Extended for NRS2 coverage |
| `F290LP_G395M/H` | 2.9 – 5.2 | 2.6 – 5.6 | Extended for NRS2 coverage |
| `CLEAR_PRISM` | 0.6 – 5.3 | 0.5 – 5.6 | Full detectable range |
| `F110W_G395H` | 0.7 – 1.3 | 0.7 – 5.6 | (Check specific science cases) |


## DJA msaexp v4 (MOS)

https://github.com/gbrammer/msaexp/blob/main/msaexp/pipeline_extended.py

* `F070LP_G140M`: [0.6, 3.3]
* `F100LP_G140M`: [0.9, 3.3]
* `F170LP_G235M`: [1.5, 5.3]
* `F290LP_G395M`: [2.6, 5.6]
* `F070LP_G140H`: [0.6, 3.3]
* `F100LP_G140H`: [0.9, 3.3]
* `F170LP_G235H`: [1.5, 5.3]
* `F290LP_G395H`: [2.6, 5.6]
* `CLEAR_PRISM`: [0.5, 5.6]


## Parlanti et al. (IFU)

https://github.com/eleonoraparlanti/nirspecIFU-extended  

* G140M/F070LP – gain only 0.90 – 0.97 µm; contamination
* G140M/F100LP – great benefit
* G235M – great benefit
* G395M – gain 5.27 –> 5.5 µm
