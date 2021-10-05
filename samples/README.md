# sample files for glafic version 2

Here are a few examples of how the code works. Note that these examples may require some input files which are also located here - so please run them right in this directory! Note that the example input files are updated (including `.dat` and `.fits` files) so they are not exactly same as those in version 1.

1. `example.input`  
 This is same as the one you can generate by `glafic -d`. Compute lensing properties and images for a compound (nfwpot+sie) system. The lenses are located at different redshifts to show its multi lens plane capability. Check output fits files to figure out the lens system. Also use `plot_point.py` and check an output PDF file to figure out how the code calculates critical curves and solves lens equations for point sources. 

2. `extend.input`  
 Fit observed arc images to derive the best-fit lens model. 

3. `point.input`
 Fitting a standard quad lens system with SIE+shear. The best-fit Hubble constant is also derived from time delay measurements.  Try both `chi2_splane 0` (image plane chi^2) and `chi2_splane 1` (source plane chi^2) to see the difference of them.

4. `clugal.input`  
 An example of cluster modeling that consists of a NFW halo and a number of member galaxies modeled by truncated isothermal ellipsoid (pseudo-Jaffe, actually). Run `plot_point.py` to check how the code resolves critical curves of small clumps.

5. `srcs.input`  
 An example of the extended source `srcs`. In the example, both lensed extended images and PSFs are produced. You can also check the effect of different PSFs and obs parameters (gain etc) by modifying parameters for `psf` in the input file.

6. `quadhost.input`  
 More complicated fitting with extended sources. This is an example of fitting quadruply imaged quasar lens system including its host galaxy component. Fitting lens, extended sources, and PSF components simultaneously.  

7. `multilens_point.input`  
 Fitting for multi lens planes with 3 SIS lenses at different redshifts. Additional point source at higher redshift. Assume that positions of all the 3 SIS lenses are observed and therefore constrained. The code tries to fit redshifts of 2 background SIS lenses as well.

8. `multilens_extend.input`  
 Similar to `multilens_point.input`, but fitting of extended sources. The foreground galaxy is also added. 

9. `point_mcmc.input`  
 Similar to `point.input`, but an example of MCMC analysis is shown. 

last modified: 10/05/2021 (based on version 2.1.2)

