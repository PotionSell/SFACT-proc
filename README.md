# SFACT-proc

Spectral processing codes for the Star Formation Across Cosmic Time (SFACT) extragalactic survey.

Processing is currently centered on WRALF (WRapped Automated Line Fitting). WRALF is, at its core, a Python wrapper for a
customized version of the automated line-fitting code ALFA (customized version: https://github.com/PotionSell/ALFA and source:
https://github.com/rwesson/ALFA).

WRALF operates on one- or two-dimensional spectral files in .fits format. It uses PyRAF (IRAF packaged as a Python module) to let the user interact with spectra via IRAF's splot tool.

To use WRALF, you will need:
1) Python 3+
2) IRAF with an existing login.cl (required by PyRAF)
3) PyRAF (http://www.stsci.edu/institute/software_hardware/pyraf)
4) the customized version of ALFA (https://github.com/PotionSell/ALFA)
5) the SFACT-proc code

To install these: 

PyRAF: use the linked instructions above.

ALFA: 
1) Create a `github` folder in your home directory
2) Move to this folder, and install the customized ALFA by cloning from its repository: 
  `>>> git clone https://github.com/PotionSell/ALFA.git`
3) Follow the installation instructions at https://github.com/PotionSell/ALFA. This involves moving to the ALFA folder, and running `make` and then `make install`

SFACT-proc:
1) Move to the `github` folder
2) Install SFACT-proc by cloning from its repository:
  `>>> git clone https://github.com/PotionSell/SFACT-proc.git`

To update either software package, simply go to its corresponding directory and use `git pull`. This will automatically clone the software from its GitHub repository.


You are now ready to process some spectra!

In standard use, the program structure is as follows:
1) the user inputs the full file path to the spectra data file
2) if available, the user inputs the full file path to the sky lines file
3) starting with the first spectrum, the user is prompted to measure one emission line with a known rest wavelength
4) this measurement is turned into a redshift fiducial to call ALFA; ALFA provides measurements of the continuum, emission line
identities, fluxes, and equivalent widths
5) the continuum and emission lines are displayed and saved via Matplotlib
6) the user can then remeasure the spectrum, or keep ALFA's estimates. If kept, the estimates and other data are packaged in 
Pandas dataframes
7) steps 1-4 are repeated for each spectrum
8) data is written into two files
9) emission line data is used to compute reddening coefficients and corrected line ratios; this is also saved to file

(more documentation soon to follow)
