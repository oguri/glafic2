#!/usr/bin/env python
import glafic
import numpy as np

glafic.version()

glafic.init(0.3, 0.7, -1.0, 0.7, 'out', -60.0, -60.0, 60.0, 60.0, 0.2, 3.0, 5, ran_seed = -10, verb = 0)

glafic.set_secondary('chi2_splane 1', verb = 0)

glafic.startup_setnum(2, 2, 1)
glafic.set_lens(1, 'nfwpot', 0.3, 7.2e14, 0.0, 0.0, 0.3, -45.0, 6.0,  0.0)
glafic.set_lens(2, 'sie',    0.5, 300.0,  2.0, 2.0, 0.2, -20.0, 0.02, 0.0)
glafic.set_extend(1, 'sersic', 1.5, 150.0, -1.0, -1.5, 0.3, 90.0, 0.8, 1.0)
glafic.set_extend(2, 'gauss',  2.0, 150.0,  1.2,  1.0, 0.2, 10.0, 0.6, 0.0)
glafic.set_point(1, 2.5, 1.0, 0.5)

glafic.setopt_lens(1, 0, 1, 0, 0, 1, 1, 0, 0)
glafic.setopt_lens(2, 0, 1, 0, 0, 1, 1, 0, 0)
glafic.setopt_extend(1, 0, 1, 1, 1, 1, 1, 0, 0)
glafic.setopt_extend(2, 0, 1, 1, 1, 1, 1, 0, 0)
glafic.setopt_point(1, 0, 1, 1)

# model_init needs to be done again whenever model parameters are changed
glafic.model_init(verb = 0)

nx_ext, ny_ext = glafic.get_nxy_ext()
print('nx_ext = %d' % nx_ext)
print('ny_ext = %d' % ny_ext)

a = glafic.calcimage(2.5, 1.0, -1.5, verb = 0)
print(a)

print('ein(id=1) = %f' % glafic.calcein_i(2.5, id = 1))
print('ein(id=1) = %f' % glafic.calcein_i(2.5, id = 2))
print('ein2 = %f' % glafic.calcein2(2.5, 0.0, 0.0))
print('kappa_ave = %f' % glafic.kappa_ave(2.5, 0.0, 0.0, 16.0))
print('kappa_cum = %f' % glafic.kappa_cum(2.5, 0.0, 0.0, 16.0))

a = glafic.point_solve(1.5, 1.0, 0.5, verb = 0)
print(a)

#glafic.readpsf('../samples/out_psf.fits', verb = 1)

# return extended image by tuple
a = glafic.writeimage()
i = 300
j = 250
print('pixel value at i = %d and j = %d is %f' % (i, j , a[j][i]))

glafic.quit()

