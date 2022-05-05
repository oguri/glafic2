#!/usr/bin/env python
import glafic
import numpy as np

glafic.version()

glafic.init_file('../samples/example.input', verb = 0)

# model_init not needed if init_file is used 
#glafic.model_init(verb = 0)

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

