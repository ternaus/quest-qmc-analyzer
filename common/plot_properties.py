__author__ = 'vladimir'

from pylab import *

w = 8
h = 8

params = {'backend': 'ps',
          'axes.labelsize': 30,
          'text.fontsize': 35,
          'legend.fontsize': 30,
          'xtick.labelsize': 30,
          'ytick.labelsize': 30,
          'figure.subplot.left': 0.10,
          'figure.subplot.bottom': 0.14,
          'figure.subplot.right': 0.97,
          'figure.subplot.top': 0.95,
          'figure.figsize': [w, h],
          'xtick.major.top': False,
}
rcParams.update(params)
ax = subplot(111)
ax.xaxis.grid(True)
ax.yaxis.grid(True)
