#!/usr/bin/env python
from __future__ import division

__author__ = 'vladimir'

'''
This script is specifically for the lattice with 3 orbitals.
y_variable:
 01 - nn correlation between 0 and 1 orbitals
 11 - nn correlation between 1 and 1 orbitals

'''
import common.get_file_list
import common.fequals
import common.divide_into_classes
import common.merge
import common.choices
import argparse
from pylab import *
import os

execfile(os.path.join(os.getcwd(), "common", "plot_properties.py"))

parser = argparse.ArgumentParser()
parser.add_argument('-t', type=float, default=1, help="hopping strength t")
parser.add_argument('-mu', type=float, help="chemical potential")
parser.add_argument('-u', type=float, help="u term")
parser.add_argument('-beta', type=float, help="inverse temperature")

parser.add_argument('-m', type=str, help="Name of the model")
parser.add_argument('-y_variable', type=str, help="variable along y axis")
parser.add_argument('-x_variable', type=str, help="variable along x axis")
parser.add_argument('-x_min', type=float, help='minimum x value')
parser.add_argument('-x_max', type=float, help='maximum x value')
parser.add_argument('-y_min', type=float, help='minimum y value')
parser.add_argument('-y_max', type=float, help='maximum y value')
parser.add_argument('-to_screen', default=False, type=bool, help='Should we print results on the screen or not.')

parser.add_argument('-legend', type=str, help='legend position. Possible values: lr, ur, ll, ul')

args = parser.parse_args(sys.argv[1:])
modelName = args.m

path = os.path.join('..', '..', 'results', modelName)

dataList = common.get_file_list.get_filelist(modelName, path)

if args.x_variable == 'rho':
  xlabel(r'$x$')
  title(r"$u = {u}$, $\beta = {beta}$".format(u=args.u, beta=args.beta), fontsize=30)
  dataList = [item for item in dataList if (common.fequals.equals(item.get_t_up(), args.t)
                                            and common.fequals.equals(item.get_u(), args.u)
                                            and common.fequals.equals(item.get_beta(), args.beta))]
  # if args.mu == 0:
  #   dataList = [item for item in dataList if (abs(item.get_rho()[0] - 1) <= item.get_rho()[1])]


#divide into classes, corresponding to different number of sites
into_nSites_dict = common.divide_into_classes.divide_into_classes(dataList, parameter='nSites')

#We choose what lattice sizes are we interested in
nSites_list = common.choices.list_choice(into_nSites_dict.keys())

for nSites in nSites_list:
  splitted = common.divide_into_classes.divide_into_classes(into_nSites_dict[nSites], parameter=args.x_variable)
  xList, yList, yErr = common.merge.merge(splitted, args.y_variable)

  xList = []
  yList = []
  yErr = []
  print into_nSites_dict[nSites]
  for temp in into_nSites_dict[nSites]:
    xList += [temp.get_rho()[0]]
    yList += [temp.get_density_correlation()[(0, 1, 0.5, 0, 0)][0]]
    yErr += [temp.get_density_correlation()[(0, 1, 0.5, 0, 0)][1]]
  title(r"$u = {u}$, $\beta = {beta}, nSites={nSites}$".format(u=args.u, beta=args.beta,
                                                               nSites=nSites), fontsize=30)
  # for key, value in temp.iteritems():
  errorbar(xList, yList, yerr=yErr, fmt='D-', label='orbits = {orbits}\n'.format(orbits='0,1'), linewidth=3)

ylabel('pairing correlation function')

if args.legend == 'lr' or args.legend == 'rl':
  legend(loc='lower right')
elif args.legend == 'ur' or args.legend == 'ru':
  legend(loc='upper right')
elif args.legend == 'll' or args.legend == 'll':
  print 'here'
  legend(loc='lower left')
elif args.legend == 'ul' or args.legend == 'lu':
  legend(loc='upper left')
else:
  legend()

try:
  ylim(ymin=args.y_min)
except:
  pass

show()
