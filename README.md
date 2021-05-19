# glafic version 2

glafic is a software for studying gravitational lensing. Use it at your own risk; the author shall not take any responsibility for loss or damage caused by the use of this software. The binary files are available at [this URL](https://www.slac.stanford.edu/~oguri/glafic/). If you use this software (or any modified version of it) for your research work, please cite the following paper:

- [M. Oguri, PASJ, 62, 1017 (2010)](https://ui.adsabs.harvard.edu/abs/2010PASJ...62.1017O/abstract) 

## Installation

glafic requires cfitsio, fftw3, and gsl to compile. Please make sure you have the library and header files in your system. If not, you can download them at [this URL for cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/), [this URL for fftw](http://www.fftw.org/index.html), and [this URL for gsl](https://www.gnu.org/software/gsl/), or can install with homebrew on macOS.

## Updates from version 1

Because of new implementation of multiple lens plane gravitational lensing, input files for version 2 are not compatible with those for version 1. Below I list important updates and tips for switching from version 1 to version 2.

* Since multiple lens redshifts are now allowed, in version 2 lens redshifts are specified as one of model parameters of each lens, rather than as a primary parameter. This means that the code uses 8 parameters to characterize each lens component, in contrast to 7 parameters adopted in version 1. Input files need to be modified to reflect this change. Since the labeling of model parameters of lens components is also changed (e.g., in version 2, 3rd and 4th model parameters refer to the central position of the lens component, rather than 2nd and 3rd model parameters used in version 1), prior files used in optimization should also be modified in most cases.

* In version 2 routines from Numerical Recipes are no longer used in order to make the source code publicly available. Several routines such as numerical integrations, random number generation, optimizations, and interpolations, are changed from version 1 to version 2. In addition, some constants such as speed of light are also updated to include more digits. Because of these changes, output values can be slightly different between version 1 and 2, even for the exactly same input parameters. The differences should be small enough in most cases though.

* Python interface is added such that some glafic routines can be called from python. The module `glafic.so` is created by running `make all`. Some examples can be found in `test_python` directory. 

## Feedback

If you have any questions, suggestions, or comments regarding glafic, please contact [Masamune Oguri](https://oguri.github.io/).

## History (only version 2)

| date       | version | comments |
|:---        |:---     |:---      |
| 2021.05.19 | 2.0.0   | first release of version 2 |
| 2021.05.03 | 2.0b6   | use a fast analytic calculation for `pow` by default |
| 2021.04.30 | 2.0b5   | `jaffe` extended source model added |
| 2021.04.27 | 2.0b2   | python interface added |
| 2021.03.23 | 2.0b1   | uploaded beta 1 version to github |
