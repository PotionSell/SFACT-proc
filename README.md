# SFACT-proc

Spectral processing codes for the Star Formation Across Cosmic Time (SFACT) extragalactic survey.

Processing is currently centered on WRALF (WRapped Automated Line Fitting). WRALF is, at its core, a Python wrapper for a
customized version of the automated line-fitting code ALFA (customized version: https://github.com/PotionSell/ALFA and source:
https://github.com/rwesson/ALFA).

WRALF operates on one- or two-dimensional spectral files in .fits format. It uses Pyraf (IRAF packaged as a Python module) to let the user interact with spectra via IRAF's splot tool.

In standard use, the program structure is as follows:
1) starting with the first spectrum, the user is prompted to measure one emission line with a known rest wavelength
2) this measurement is turned into a redshift fiducial to call ALFA; ALFA provides measurements of the continuum, emission line
identities, fluxes, and equivalent widths
3) the continuum and emission lines are displayed and saved via Matplotlib
4) the user can then remeasure the spectrum, or keep ALFA's estimates. If kept, the estimates and other data are packaged in 
Pandas dataframes
5) steps 1-4 are repeated for each spectrum
6) data is written into two files
7) emission line data is used to compute reddening coefficients and corrected line ratios; this is also saved to file

(more documentation soon to follow)
