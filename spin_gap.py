# !/usr/bin/env python
from __future__ import division

__author__ = 'vladimir'

'''
This script will be used to plot the spin gap versus rho

min_beta - minimum beta value that is used. rest will be averaged.

Uses file where defined minimum beta that is used for extrapolation for different rho
'''
import common.get_file_list
import common.fequals
import common.divide_into_classes
import common.merge
import common.choices
import common.extrapolation

import argparse
from pylab import *
import os

execfile(os.path.join(os.getcwd(), "common", "plot_properties.py"))

parser = argparse.ArgumentParser()
parser.add_argument('-t', type=float, default=1, help="hopping strength t")
parser.add_argument('-u', type=float, help="u term")

parser.add_argument('-file', type=str, help="file where minimum beta for different rho is stored")

parser.add_argument('-m', type=str, help="Name of the model")
parser.add_argument('-x_min', type=float, help='minimum x value')
parser.add_argument('-x_max', type=float, help='maximum x value')
parser.add_argument('-y_min', type=float, help='minimum y value')
parser.add_argument('-y_max', type=float, help='maximum y value')
parser.add_argument('-to_screen', default=False, type=bool, help='Should we print results on the screen or not.')

parser.add_argument('-legend', type=str, help='legend position. Possible values: lr, ur, ll, ul')

args = parser.parse_args(sys.argv[1:])
modelName = args.m

execfile('settings.py')
execfile(args.file)

path = os.path.join(folder_with_different_models, modelName)

dataList = common.get_file_list.get_filelist(modelName, path)

#Clean up datafiles that are for sure bad. No moves were performed.
dataList = [item for item in dataList if
            ((item.get_global_sites() > 0 and item.get_global_accept > 0) or item.get_global_sites() == 0)]

# Clean up datafiles that are for sure bad. Electron density rho is bounded by 0 and 2.
dataList = [item for item in dataList if (item.get_rho()[0] < 2.1)]

#Not using files, that they have mu = 0 and 0 global moves, when u is not zero.
# dataList = [item for item in dataList if
#             ((
#              item.get_u() != 0 and item.get_mu_up() == 0 and item.get_global_sites() > 0) or item.get_u() == 0 or item.get_mu_up() != 0)]

if args.m == 'Lieb' or args.m == 'square':
  bipartite = True
else:
  bipartite = False

# if bipartite:
#   #We know from symmetry of the bipartite lattice that at mu = 0, rho = 1, I will assume that 0.97 - 1.03 is ok
#   if args.rho == 1:
#     dataList = [item for item in dataList if (item.get_mu_up() == 0)]
#     dataList = [item for item in dataList if (abs(item.get_rho()[0] - 1) <= 0.03)]

# We need to remove datapoints if s-wave errorbars > 20%
dataList = [item for item in dataList if (abs(item.get_s_wave()[1] / item.get_s_wave()[0]) < 0.2)]


# title(r"{modelname}, $\rho = {rho}$, $u = {u}$, $min \beta = {beta}$, $t={t}$".format(u=args.u, rho=args.rho, beta=args.min_beta,
#                                                                            modelname=args.m, t=args.t), fontsize=30)

#filter files with respect to U and t
dataList = [item for item in dataList if (common.fequals.equals(item.get_t_up(), args.t)
                                          and common.fequals.equals(item.get_u(), args.u))]
# and (abs(item.get_rho()[0] - args.rho)) < 0.02 + item.get_rho()[1]
# and (item.get_beta() >= args.min_beta))]

#Now we go thorough the rhoList, divide into the classes with respect to the rhoList and keep only beta that we are interested in

cList = []
cErr = []

for rho in sorted(rhoDict.keys()):
  print rho
  rho_data = [item for item in dataList if ((abs(item.get_rho()[0] - rho)) < 0.02 + item.get_rho()[1]
                                            and item.get_beta() >= rhoDict[rho])]
  #divide into classes, corresponding to different number of sites
  into_nSites_dict = common.divide_into_classes.divide_into_classes(rho_data, parameter='nSites')

  #Each lattice will give one value

  xList = []
  yList = []
  yErr = []

  xList, yList, yErr = common.merge.merge(into_nSites_dict, 's-wave')

  #renormalization
  for i in range(len(xList)):
    yList[i] = yList[i] / xList[i]
    yErr[i] = yErr[i] / yErr[i]
    xList[i] = 1 / math.sqrt(xList[i])

  print xList
  print yList
  print yErr
  approxx, approxy, m, c, std_err = common.extrapolation.approx(xList, yList, yErr)

  c = max(c, 0)
  cList += [math.sqrt(c)]
  cErr += [std_err]

errorbar(sorted(rhoDict.keys()), cList, yerr=cErr, fmt='D-', linewidth=3,
         markersize=15)
ylabel(r'$\Delta_0$')
xlabel(r'$\rho$')

if args.to_screen:
  print 'xList = ', sorted(rhoDict.keys())
  print 'yList = ', cList
  print 'yErr = ', cErr

if '-x_min' in sys.argv:
  xlim(xmin=args.x_min)

show()
