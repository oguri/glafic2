# glafic version 2

glafic is a software for studying gravitational lensing. Use it at your own risk; the author shall not take any responsibility for loss or damage caused by the use of this software. The binary files are available at [this URL](https://drive.google.com/drive/folders/1AJwbOiHbaV65_cEo3qbz3c4s3tvD7rS8?usp=share_link). If you use this software (or any modified version of it) for your research work, please cite the following paper:

- [M. Oguri, PASJ, 62, 1017 (2010)](https://ui.adsabs.harvard.edu/abs/2010PASJ...62.1017O/abstract) 

and if `anfw` or `ahern` model is used it would be appropriate to also cite the following paper:

- [M. Oguri, PASP, 133, 074504 (2021)](https://ui.adsabs.harvard.edu/abs/2021PASP..133g4504O/abstract) 

## Installation

glafic requires cfitsio, fftw3, and gsl to compile. Please make sure you have the library and header files in your system. If not, you can download them at [this URL for cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/), [this URL for fftw](http://www.fftw.org/index.html), and [this URL for gsl](https://www.gnu.org/software/gsl/), or can install with homebrew on macOS. After installing those libraries, you can simply run
```
make
```
to compile glafic and create an executable file. 

## Updates from version 1

Because of new implementation of multiple lens plane gravitational lensing, input files for version 2 are not compatible with those for version 1. Below I list important updates and tips for switching from version 1 to version 2.

* Since multiple lens redshifts are now allowed, in version 2 lens redshifts are specified as one of model parameters of each lens, rather than as a primary parameter. This means that the code uses 8 parameters to characterize each lens component, in contrast to 7 parameters adopted in version 1. Input files need to be modified to reflect this change. Since the labeling of model parameters of lens components is also changed (e.g., in version 2, 3rd and 4th model parameters refer to the central position of the lens component, rather than 2nd and 3rd model parameters used in version 1), prior files used in optimization should also be modified in most cases.

* In version 2 routines from Numerical Recipes are no longer used in order to make the source code publicly available. Several routines such as numerical integrations, random number generation, optimizations, and interpolations, are changed from version 1 to version 2. In addition, some constants such as speed of light are also updated to include more digits. Because of these changes, output values can be slightly different between version 1 and 2, even for the exactly same input parameters. The differences should be small enough in most cases though.

* Python interface is added such that some glafic routines can be called from python (see below).

## Python interface

glafic can be run with python. Some examples can be found in `test_python` directory. To run glafic on python, please first try
```
make python
```
in your python environment to see if a library file is successfully created as `python/glafic/glafic.so`. Then you can run
```
pip install .
```
in the directory where `setup.py` is located to install glafic. Once glafic is installed, you can use glafic in your python scripts simply by importing the module as `import glafic`. 

## Feedback

If you have any questions, suggestions, or comments regarding glafic, please contact [Masamune Oguri](https://oguri.github.io/).

## History (only version 2)

| date       | version | comments |
|:---        |:---     |:---      |
| 2024.10.02 | 2.1.10  |  another option for `chi2_usemag` added |
| 2024.05.27 | 2.1.9   |  `gaupot` lens model added |
| 2023.12.06 | 2.1.8   | `getpar_*` commands added for the python interface |
| 2023.11.30 | 2.1.7   | minor bug fix |
| 2023.05.04 | 2.1.6   | `optpoint` and `optextend` added for the python interface |
| 2022.10.29 | 2.1.5   | `coord2xy` and `xy2coord` added, bug fix |
| 2022.09.26 | 2.1.4   | `crline` lens model added, `flag_outpot` added |
| 2022.05.05 | 2.1.3   | `init_file` added for the python interface |
| 2021.08.31 | 2.1.2   | `flatfix` added, bug fix |
| 2021.07.27 | 2.1.1   | `mcmc_kapcum` and `mcmc_findimg` command added |
| 2021.06.22 | 2.1.0   | `anfw` and `ahern` lens models added, `c2calc` and `reset_obs_point` command added |
| 2021.05.19 | 2.0.0   | first release of version 2 |
| 2021.05.03 | 2.0b6   | use a fast analytic calculation for `pow` by default |
| 2021.04.30 | 2.0b5   | `jaffe` extended source model added |
| 2021.04.27 | 2.0b2   | python interface added |
| 2021.03.23 | 2.0b1   | uploaded beta 1 version to github |
