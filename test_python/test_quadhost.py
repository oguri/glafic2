#!/usr/bin/env python
import glafic
import numpy as np

glafic.init(0.3, 0.7, -1.0, 0.7, 'out', -5.0, -5.0, 5.0, 5.0, 0.1, 0.5, 6, verb = 0)

glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('obs_gain       3.0', verb = 0)
glafic.set_secondary('obs_ncomb      10', verb = 0)
glafic.set_secondary('psfconv_size   3.0', verb = 0)
glafic.set_secondary('flag_extnorm   1', verb = 0)

glafic.startup_setnum(2, 3, 0)
glafic.set_lens(1, 'sie', 0.3, 315.0, 0.0, 0.0, 0.3,  37.0, 0.0, 0.0)
glafic.set_lens(2, 'pert', 0.3,  2.0, 0.0, 0.0, 0.042, 60.0, 0.0, 0.0)
glafic.set_extend(1,'sersic', 0.3, 4.9e5,  0.0,  0.0, 0.3,   37.0, 1.0, 4.0)
glafic.set_extend(2, 'sersic', 2.0, 1.2e4,  0.21, -0.09, 0.3, -10.0, 0.5, 2.0)
glafic.set_extend(3, 'point',  2.0, 2.1e4,  0.21, -0.09, 0.0,   0.0, 0.0, 0.0)
glafic.set_psf(0.22, 0.07, 15.0, 3.0, 0.8, 0.04, -20.0, 2.5, 0.35)

#glafic.setopt_lens(1, 0, 1, 0, 0, 0, 0, 0, 0)
#glafic.setopt_lens(2, 0, 0, 0, 0, 1, 1, 0, 0) 
#glafic.setopt_extend(1, 0, 1, 0, 0, 0, 0, 1, 0)
#glafic.setopt_extend(2, 0, 1, 1, 1, 0, 0, 1, 0) 
#glafic.setopt_extend(3, 0, 1, 1, 1, 0, 0, 0, 0) 
#glafic.setopt_psf(1, 0, 0, 0, 1, 0, 0, 0, 1) 

# model_init needs to be done again whenever model parameters are changed
glafic.model_init(verb = 0)

glafic.readobs_extend('../samples/obs_quadhost.fits', verb = 0 )
#glafic.optimize();
print('chi2 = %e' % glafic.c2calc())
  
glafic.set_lens(1, 'sie', 0.3, 320.0, 0.0, 0.0, 0.3,  37.0, 0.0, 0.0)
glafic.model_init(verb = 0)

print('chi2 = %e' % glafic.c2calc())

glafic.quit()

