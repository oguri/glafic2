# glafic version 2

glafic is a software for studying gravitational lensing. Use it at your own risk; the author shall not take any responsibility for loss or damage caused by the use of this software. The binary files are available at [this URL](https://www.slac.stanford.edu/~oguri/glafic/).

## Installation

glafic requires cfitsio, fftw3, and gsl to compile. Please make sure you have the library and header files in your system. If not, you can download them at [this URL for cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/), [this URL for fftw](http://www.fftw.org/index.html), and [this URL for gsl](https://www.gnu.org/software/gsl/).

## Updates from version 1

Because of new implementation of multiple lens plane gravitational lensing, input files for version 2 are not compatible with those for version 1. Below I list important updates and tips for switching from version 1 to version 2.

* Since multiple lens redshifts are now allowed, in version 2 lens redshifts are specified as one of model parameters of each lens, rather than as a primary parameter. This means that the code uses 8 parameters to characterize each lens component, in contrast to 7 parameters adopted in version 1. Input files need to be modified to reflect this change. Since the labeling of model parameters of lens components is also changed (e.g., in version 2, 3rd and 4th model parameters refer to the central position of the lens component, rather than 2nd and 3rd model parameters used in version 1), prior files used in optimization should also be modified in most cases.

* In version 2 routines from Numerical Recipes are no longer used in order to make the source code publicly available. Several routines such as numerical integrations, random number generation, optimizations, and interpolations, are changed from version 1 to version 2. In addition, some constants such as speed of light are also updated to include more digits. Because of these changes, output values can be slightly different between version 1 and 2, even for the exactly same input parameters. The differences should be small enough in most cases though.

## Feedback

If you have any questions, suggestions, or comments regarding glafic, please contact [Masamune Oguri](https://oguri.github.io/).

## History (only version 2)

| date       | version | comments |
|:---        |:---     |:---      |
| 2021.03.23 | 2.0b1   | uploaded beta 1 version to github |
