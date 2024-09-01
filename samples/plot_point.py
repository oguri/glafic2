#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.ticker import StrMethodFormatter, NullFormatter

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.size'] = 20
mpl.rcParams['legend.handlelength'] = 3.0
mpl.rcParams['axes.linewidth'] = 1.0
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['xtick.minor.width'] = 1.0
mpl.rcParams['xtick.major.size'] = 5.0
mpl.rcParams['xtick.minor.size'] = 3.0
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['ytick.major.width'] = 1.0
mpl.rcParams['ytick.minor.width'] = 1.0
mpl.rcParams['ytick.major.size'] = 5.0
mpl.rcParams['ytick.minor.size'] = 3.0
mpl.rcParams['ytick.right'] = True

plt.rcParams['ps.useafm'] = True
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = '\\usepackage{sfmath}'

rfile = 'out_crit.dat'
data = np.loadtxt(rfile, comments = '#')

ix1 = data[: ,0]
iy1 = data[: ,1]
sx1 = data[: ,2]
sy1 = data[: ,3]
ix2 = data[: ,4]
iy2 = data[: ,5]
sx2 = data[: ,6]
sy2 = data[: ,7]

rfile = 'out_mesh.dat'
data = np.loadtxt(rfile, comments = '#')

ix3 = data[: ,0]
iy3 = data[: ,1]
sx3 = data[: ,2]
sy3 = data[: ,3]
ix4 = data[: ,4]
iy4 = data[: ,5]
sx4 = data[: ,6]
sy4 = data[: ,7]

rfile = 'out_point.dat'
data = np.loadtxt(rfile, comments = '#')

nimg = 0
sxx = []
syy = []
ixx = []
iyy = []
for i in range(data.shape[0]):
    if nimg == 0:
        nimg = data[i, 0]
        sxx.append(data[i, 2])
        syy.append(data[i, 3])
    else:
        nimg -= 1
        ixx.append(data[i ,0])
        iyy.append(data[i ,1])

plt.figure(figsize = (8, 16))

plt.subplot(2, 1, 1)

xmin = -30
xmax = 30
ymin = -30
ymax = 30

plt.plot([ix3, ix4], [iy3, iy4], '-', color = 'grey', zorder = 0, lw = 0.5)
plt.plot([ix1, ix2], [iy1, iy2], '-', color = 'blue', zorder = 1)

plt.scatter(ixx, iyy, s = 30, marker = '^', color = 'red', zorder = 2)

ax = plt.gca()
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])

plt.title('image plane') 
plt.xlabel('$\\theta_1$ [arcsec]')
plt.ylabel('$\\theta_2$ [arcsec]')
plt.subplots_adjust(top = 0.95)
plt.subplots_adjust(bottom = 0.15)
plt.subplots_adjust(left = 0.15)

plt.subplot(2, 1, 2)

xmin = -15
xmax = 15
ymin = -15
ymax = 15

plt.plot([sx3, sx4], [sy3, sy4], '-', color = 'grey', zorder = 0, lw = 0.5)
plt.plot([sx1, sx2], [sy1, sy2], '-', color = 'blue', zorder = 1)

plt.scatter(sxx, syy, s = 30, marker = '^', color = 'red', zorder = 2)

ax = plt.gca()
plt.xlim([xmin, xmax])
plt.ylim([ymin, ymax])

plt.title('source plane') 
plt.xlabel('$\\beta_1$ [arcsec]')
plt.ylabel('$\\beta_2$ [arcsec]')
plt.subplots_adjust(top = 0.95)
plt.subplots_adjust(bottom = 0.15)
plt.subplots_adjust(left = 0.15)

ofile = 'plot_point'

#plt.savefig(ofile + '.png', dpi = 150, transparent = True)
plt.savefig(ofile + '.pdf', bbox_inches = 'tight')
