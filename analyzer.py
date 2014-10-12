#!/usr/bin/env python
from __future__ import division

__author__ = 'vladimir'

'''
This script will be used to analyze the data generated by the QUEST simulation

x_variables:
  u - fixed mu, beta
  mu - fixed u, beta
  beta - fixed u, mu
  T - fixed u, mu
  x - fixed u, mu, beta
  1L - plot versus 1 / L. fixed u, mu, beta
  rho - plot versus rho. fixed u, beta
  num_sites - plot versus number of sites. fixed rho, u, beta

y_variables:
  energy - total energy
  energy_hop - hopping energy
  m2 - square of the magnetisation
  01 - density-density correlation between 0 and 1 orbital within unit cell
  11 - density-density correlation between 0 and 1 orbital within unit cell
  m0_squared - square of the magnetisation on the 0 orbital
  m1_squares - square of the magnetisation on the 1 orbital
  C - specific heat
  s-wave_rescaled - L^(-7/4) * Ps
  sign_up_down - sign_up times sign_down
'''
import common.get_file_list
import common.fequals
import common.divide_into_classes
import common.merge
import common.choices
from pylab import *
import argparse
import os
import time

start_time = time.time()

execfile(os.path.join(os.getcwd(), "common", "plot_properties.py"))

parser = argparse.ArgumentParser()
parser.add_argument('-t', type=float, default=1, help="hopping strength t")
parser.add_argument('-t1', type=float, default=0, help="extra hopping strength t1")
parser.add_argument('-mu', type=float, help="chemical potential")
parser.add_argument('-u', type=float, help="u term")
parser.add_argument('-beta', type=float, help="inverse temperature")

parser.add_argument('-m', type=str, help="Name of the model")
parser.add_argument('-y_variable', type=str, help="""variable along y axis can be
sign - average sign of the determinant
sign_up - average sign of the determinant for spin up electrons
sign_down - average sign of the determinant for spin down electrons
m2_rho - m2 divided by rho.
sign_up_down - square sign_up * sign_down""")

parser.add_argument('-x_variable', type=str, help="variable along x axis")
parser.add_argument('-x_min', type=float, help='minimum x value')
parser.add_argument('-x_max', type=float, help='maximum x value')
parser.add_argument('-y_min', type=float, help='minimum y value')
parser.add_argument('-y_max', type=float, help='maximum y value')
parser.add_argument('-to_screen', default=False, type=bool, help='Should we print results on the screen or not.')

parser.add_argument('-T', action="store_true",
                    help="if we need graph versus temperature. By default it is versus inverse temperature.")
parser.add_argument('-legend', type=str, help='legend position. Possible values: lr, ur, ll, ul')

parser.add_argument('-filter', type=bool, default=True, help='Do we filter the data or not. Default True')
parser.add_argument('-vline', type=float, help='Adds vertical line to the plot at a given position')
parser.add_argument('-hline', type=float, help='Adds horizontal line to the plot at a given position')

args = parser.parse_args(sys.argv[1:])
modelName = args.m

execfile('settings.py')

path = os.path.join(folder_with_different_models, modelName)

if args.m == 'chain':
  dimension = 1
else:
  dimension = 2

dataList = common.get_file_list.get_filelist(modelName, path, dimension=dimension)

#Clean up datafiles that are for sure bad. No moves were performed.
dataList = (item for item in dataList if
                ((item.get_global_sites() > 0 and item.get_global_accept > 0) or item.get_global_sites() == 0))

# Clean up datafiles that are for sure bad. Electron density rho is bounded by 0 and 2.


dataList = (item for item in dataList if (0 <= item.get_rho()[0] <= 2))

# Clean up datafiles that have average sign less than 0.05

# dataList = (item for item in dataList if (item.get_sign()[0] > 0.05))

#Not using files, that they have mu = 0 and 0 global moves, when u is not zero.
# dataList = [item for item in dataList if
#             ((
#              item.get_u() != 0 and item.get_mu_up() == 0 and item.get_global_sites() > 0) or item.get_u() == 0 or item.get_mu_up() != 0)]

bipartite_list = ['Lieb', 'square', 'honeycomb', 'chain', '9_16_depleted']
if args.m in bipartite_list:
  bipartite = True
else:
  bipartite = False

# if bipartite:
#   #We know from symmetry of the bipartite lattice that at mu = 0, rho = 1, I will assume that 0.97 - 1.03 is ok
#   if args.rho == 1:
#     dataList = [item for item in dataList if (item.get_mu_up() == 0)]
#     dataList = [item for item in dataList if (abs(item.get_rho()[0] - 1) <= 0.03)]

if args.filter:
  # We need to remove datapoints if s-wave errorbars > 20%
  dataList = (item for item in dataList if (abs(item.get_s_wave()[1] / item.get_s_wave()[0]) < 0.2))

# Filter t
dataList = (item for item in dataList if (common.fequals.equals(item.get_t_up()[0], args.t)))
#Filter t'
if 'anisotropic' in args.m:
  dataList_temp = []
  for item in dataList:
	  if len(item.get_t_up()) == 1:
		  assert args.t1 == args.t
		  dataList_temp += [item]
  dataList = dataList_temp


  # dataList = (item for item in dataList if (common.fequals.equals(item.get_t_up()[1], args.t1)))

if args.x_variable == 'u':
  xlabel(r'$U[t]$')
  title(r"{modelname}, $\rho = {rho}$, $\beta = {beta}$, $t={t}$".format(beta=args.beta, rho=args.rho, modelname=args.m, t=args.t),
        fontsize=30)
  dataList = (item for item in dataList if (common.fequals.equals(item.get_beta(), args.beta)
                                            and common.fequals.equals(item.get_mu_up(), args.mu)))
elif args.x_variable == 'mu':
  xlabel(r'$\mu[t]$')
  title(r"{modelname}, $U = {u}$, $\beta = {beta}$, $t={t}$".format(beta=args.beta, u=args.u, modelname=args.m, t=args.t), fontsize=30)
  dataList = (item for item in dataList if
              (common.fequals.equals(item.get_u(), args.u) and common.fequals.equals(item.get_beta(), args.beta)))

elif args.x_variable == 'rho':
  xlabel(r'$\rho$')
  title(r"{modelname}, $U = {u}$, $\beta = {beta}$, $t={t}$".format(beta=args.beta, u=args.u, modelname=args.m,
                                                                    t=[args.t, args.t1]), fontsize=30)
  dataList = (item for item in dataList if (common.fequals.equals(item.get_u(), args.u)
                                            and common.fequals.equals(item.get_beta(), args.beta)))
  if args.m == 'Lieb' or args.m == 'kagome_anisotropic':
    axvline(x=2 / 3, linestyle='--', color='black', linewidth=3)
    axvline(x=4 / 3, linestyle='--', color='black', linewidth=3)

elif args.x_variable == 'beta':
  xlabel(r'$\beta[t]$')
  title(r"{modelname}, $\mu = {mu}$, $U = {u}$, $t={t}$".format(u=args.u, mu=args.mu, modelname=args.m, t=args.t),
        fontsize=30)
  dataList = (item for item in dataList if (common.fequals.equals(item.get_u(), args.u)
                                            and (common.fequals.equals(item.get_mu_up(), args.mu))))

elif args.x_variable == 'T':
  xlabel(r'$T[t]$')
  title(r"{modelname}, $mu = {mu}$, $U = {u}$, $t={t}$".format(u=args.u, mu=args.mu, modelname=args.m, t=args.t),
        fontsize=30)
  dataList = [item for item in dataList if (common.fequals.equals(item.get_u(), args.u)
                                            and common.fequals.equals(item.get_mu_up(), args.mu))]
elif args.x_variable == '1L':
  xlabel(r'$1 / L$')
  title(r"{modelname}, $\rho = {rho}$, $u = {u}$, $\beta = {beta}$, $t={t}$".format(u=args.u, rho=args.rho, beta=args.beta,
                                                                           modelname=args.m, t=args.t), fontsize=30)
  dataList = [item for item in dataList if (common.fequals.equals(item.get_u(), args.u)
                                            and common.fequals.equals(item.get_rho()[0], args.rho)
                                            and common.fequals.equals(item.get_beta(), args.beta))]

elif args.x_variable == 'num_sites':
  xlabel(r'number of sites')
  title(
    r"{modelname}, $\rho = {rho}$, $u = {u}$, $\beta = {beta}$, $t={t}$".format(u=args.u, rho=args.rho, beta=args.beta,
                                                                                modelname=args.m, t=args.t),
    fontsize=30)
  dataList = [item for item in dataList if (common.fequals.equals(item.get_u(), args.u)
                                            and common.fequals.equals(item.get_rho()[0], args.rho)
                                            and common.fequals.equals(item.get_beta(), args.beta))]


#divide into classes, corresponding to different number of sites
into_nSites_dict = common.divide_into_classes.divide_into_classes(dataList, parameter='shape')

print 'waste of time = ', time.time() - start_time

#We choose what lattice sizes are we interested in
shape_list = common.choices.list_choice(into_nSites_dict.keys())

for shape in shape_list:
  splitted = common.divide_into_classes.divide_into_classes(into_nSites_dict[shape], parameter=args.x_variable)
  xList, yList, yErr = common.merge.merge(splitted, args.y_variable)

  if args.x_variable == 'T':
    xList = [1 / tx for tx in xList]

  errorbar(xList, yList, yerr=yErr, fmt='D-',
           label=r'${N} = {nx} \times {ny}$'.format(N=shape[2], nx=shape[0], ny=shape[1]), linewidth=3,
           markersize=15)
  if args.to_screen:
    print
    print 'N = ', shape[2]
    print 'nx = ', shape[0]
    print 'ny = ', shape[1]
    print 'xList = ', xList
    print 'yList = ', yList
    print 'yErr = ', yErr

if args.y_variable == 'energy':
  ylabel(r'$Energy$')

if args.y_variable == 'energy_hop':
  ylabel(r'$Energy_{hop}$')

if args.y_variable == 'X_F':
  ylabel(r'$struct_xx_f$')

if args.y_variable == 'rho':
  ylabel(r'$\rho$')

if args.y_variable == 'DO':
  ylabel(r'$\left<N_{up} N_{down}\right>$')

if args.y_variable == 's-wave':
  ylabel(r'$P_s$')

if args.y_variable == 's-wave_rescaled':
  ylabel(r'$L^{-7/4} P_s$')

if args.y_variable == '01':
  ylabel(r'$\left<n_0 n_1\right>[t]$')

if args.y_variable == '12':
  ylabel(r'$\left<n_1 n_2\right>[t]$')

if args.y_variable == '00':
  ylabel(r'$\left<n_{0,(0,0)} n_{0, (1, 0)}\right>[t]$')

if args.y_variable == 'm2':
  ylabel(r'$\left<m^2\right>$')

if args.y_variable == 'm2_rho':
  ylabel(r'$\left<m^2 / \rho\right>$')

if args.y_variable == 'sign_normalized':
  ylabel(
    r'$\left<sign_{\uparrow} sign_{\downarrow} \right> - \left< sign_{\uparrow} \right> \left< sign_{\downarrow} \right>$')

if args.y_variable == 'm0_squared':
	ylabel(r'$\left<m_0^2 \right>[$')

if args.y_variable == 'm1_squared':
	ylabel(r'$\left<m_1^2 \right>$')

if args.y_variable == 'sign':
    ylabel(r'$\left< {\rm sign} \right>$')

if args.y_variable == 'sign_up':
  ylabel(r'$\left<sign_{\uparrow} \right>$')


elif args.y_variable == 'C':
  ylabel(r'$C$')

elif args.y_variable == 'sign_up_down':
    ylabel(r'$\left< {\rm sign}_{\uparrow}\right> \left< S_{\downarrow} \right>$')


if args.legend == 'lr' or args.legend == 'rl':
  legend(loc='lower right', fancybox=True, shadow=True)
elif args.legend == 'ur' or args.legend == 'ru':
  legend(loc='upper right', fancybox=True, shadow=True)
elif args.legend == 'll' or args.legend == 'll':
  print 'here'
  legend(loc='lower left', fancybox=True, shadow=True)
elif args.legend == 'ul' or args.legend == 'lu':
  legend(loc='upper left', fancybox=True, shadow=True)
elif args.legend == 'best':
  legend(loc='best', fancybox=True, shadow=True)
else:
  legend(fancybox=True, shadow=True)

if '-x_min' in sys.argv:
  xlim(xmin=args.x_min)

if '-vline' in sys.argv:
    axvline(x=args.vline, linestyle='-', color='black', linewidth=2)
if '-hline' in sys.argv:
    axhline(y=args.hline, linestyle='-', color='black', linewidth=2)
show()
