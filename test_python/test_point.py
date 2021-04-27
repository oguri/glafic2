#!/usr/bin/env python
import glafic
import numpy as np

glafic.init(0.3, 0.7, -1.0, 0.7, 'out', -5.0, -5.0, 5.0, 5.0, 0.02, 0.5, 5, verb = 0)

glafic.set_secondary('chi2_splane 0', verb = 0)
glafic.set_secondary('chi2_checknimg 1', verb = 0)
glafic.set_secondary('chi2_restart   -1', verb = 0)
glafic.set_secondary('chi2_usemag    0', verb = 0)
glafic.set_secondary('hvary          1', verb = 0)

glafic.startup_setnum(2, 0, 1)
glafic.set_lens(1, 'sie', 0.5, 300.0,  0.0, 0.0, 0.35,   0.0, 0.0, 0.0)
glafic.set_lens(2, 'pert', 0.5,   2.0,  0.0, 0.0, 0.05,  60.0, 0.0, 0.0)
glafic.set_point(1, 2.0, -0.15, 0.05)

glafic.setopt_lens(1, 0, 1, 1, 1, 1, 1, 0, 0)
glafic.setopt_lens(2, 0, 0, 0, 0, 1, 1, 0, 0)
glafic.setopt_point(1, 0, 1, 1)

# model_init needs to be done again whenever model parameters are changed
glafic.model_init(verb = 0)

glafic.readobs_point('../samples/obs_point.dat')
glafic.parprior('../samples/prior_point.dat')
#glafic.optimize()
print('chi2 = %e' % glafic.c2calc())

glafic.quit()

