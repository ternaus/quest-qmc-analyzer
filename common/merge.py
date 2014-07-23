from __future__ import division

__author__ = 'Vladimir Iglovikov'

import numpy
import math


def merge(dict_to_merge, variable):
  '''
  @param dict_to_merge:
    dictionary with key - some parameter, values - list of measurements
    that should be merged into
  @param variable: variable that we need to merge
    Energy
    struct_xx_f - ferromagnetic structure factor
    rho
    DO
    s-wave
    Energy_hop
    m2
    bc - binder cumulant
    01 - density-density correlation function between 0 and 1 orbitals within unit cell
    12 - density-density correlation function between 1 and 2 orbitals within unit cell
    01_normalized - density-density correlation function between 0 and 1 orbitals within unit cell - <n0> <n1>
    12_normalized - density-density correlation function between 1 and 2 orbitals within unit cell - <n1> <n2>
    00 - density-density correlation function between orbital 0 in the (0,0) unit cell and orbital 0 in the (1,0) unit cell

    n0 - density on the 0 orbital
    n1 - density on the 1 orbital
    m0_squared - square of the magnetisation on the 0 orbital
    m1_squared - square of the magnetisation on the 1 orbital
    C - specific heat
    m2_rho - magnetisation squared divided by rho
  @return:
    xList - sorted list of the keys
    yList - list of corresponding y values
    yErr - list of corresponding error bars
  '''

  xList = dict_to_merge.keys()
  xList.sort()
  yList = []
  yErr = []

  if variable == 'Energy':
    for x in xList:
      yList += [numpy.mean([item.get_energy()[0] for item in dict_to_merge[x]])]
      temp = [item.get_energy()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'struct_xx_f':
    for x in xList:
      yList += [numpy.mean([item.get_struct_XX_F()[0] for item in dict_to_merge[x]])]
      temp = [item.get_struct_XX_F()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'rho':
    for x in xList:
      rho_average = numpy.mean([item.get_rho()[0] for item in dict_to_merge[x]])
      yList += [rho_average]
      temp = [item.get_rho()[1] ** 2 for item in dict_to_merge[x]]
      if sum(temp) == 0:
        temp = [(item.get_rho()[0] - rho_average) ** 2
                for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'DO':
    for x in xList:
      DO = numpy.mean([item.get_DO()[0] for item in dict_to_merge[x]])
      yList += [DO]
      temp = [item.get_DO()[1] ** 2 for item in dict_to_merge[x]]
      if sum(temp) == 0:
        temp = [(item.get_DO()[0] - DO) ** 2
                for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 's-wave':
    for x in xList:
      yList += [numpy.mean([item.get_s_wave()[0] for item in dict_to_merge[x]])]
      temp = [item.get_s_wave()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 's-wave_rescaled':
    for x in xList:
      yList += [numpy.mean([math.pow(item.get_L(), -7.0 / 4.0) * item.get_s_wave()[0] for item in dict_to_merge[x]])]
      temp = [(math.pow(item.get_L(), -7.0 / 4.0) * item.get_s_wave()[1]) ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'Energy_hop':
    for x in xList:
      E_hop_average = numpy.mean([item.get_energy_hop()[0] for item in dict_to_merge[x]])
      yList += [E_hop_average]
      temp = [item.get_energy_hop()[1] ** 2 for item in dict_to_merge[x]]
      if sum(temp) == 0:
        temp = [(item.get_energy_hop()[0] - E_hop_average) ** 2
                for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'm2':
    for x in xList:
      m2_average = numpy.mean([item.get_m2()[0] for item in dict_to_merge[x]])
      yList += [m2_average]
      temp = [item.get_m2()[1] ** 2 for item in dict_to_merge[x]]
      if sum(temp) == 0:
        temp = [(item.get_m2()[0] - m2_average) ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'm2_rho':
    for x in xList:
      m2_rho_average = numpy.mean([item.get_m2()[0] / item.get_rho()[0] for item in dict_to_merge[x]])
      yList += [m2_rho_average]
      temp = [item.get_m2()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]


  elif variable == 'm':
    for x in xList:
      yList += [numpy.mean([item.get_m()[0] for item in dict_to_merge[x]])]
      temp = [item.get_m()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'bc':
    for x in xList:
      yList += [numpy.mean([item.get_bc()[0] for item in dict_to_merge[x]])]
      temp = [item.get_bc()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == '01':
    for x in xList:
      yList += [numpy.mean([item.get_density_correlation()[(0, 1, 0.5, 0, 0)][0] for item in dict_to_merge[x]])]
      temp = [item.get_density_correlation()[(0, 1, 0.5, 0, 0)][1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == '12':
    for x in xList:
      yList += [numpy.mean([item.get_density_correlation()[(1, 2, -0.5, 0.5, 0)][0] for item in dict_to_merge[x]])]
      temp = [item.get_density_correlation()[(1, 2, -0.5, 0.5, 0)][1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == '01_normalized':
    for x in xList:
      yList += [numpy.mean([item.get_n01_normalized()[0] for item in dict_to_merge[x]])]
      temp = [item.get_n01_normalized()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == '12_normalized':
    for x in xList:
      yList += [numpy.mean([item.get_n12_normalized()[0] for item in dict_to_merge[x]])]
      temp = [item.get_n12_normalized()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == '00':
    for x in xList:
      yList += [numpy.mean([item.get_density_correlation()[(0, 0, 1, 0, 0)][0] for item in dict_to_merge[x]])]
      temp = [item.get_density_correlation()[(0, 0, 1, 0, 0)][1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'n0':
    for x in xList:
      yList += [numpy.mean([item.get_n0()[0] for item in dict_to_merge[x]])]
      temp = [item.get_n0()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'n1':
    for x in xList:
      yList += [numpy.mean([item.get_n1()[0] for item in dict_to_merge[x]])]
      temp = [item.get_n1()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]


  elif variable == 'm0_squared':
    for x in xList:
      yList += [numpy.mean([item.get_m0_squared()[0] for item in dict_to_merge[x]])]
      temp = [item.get_m0_squared()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'm1_squared':
    for x in xList:
      yList += [numpy.mean([item.get_m1_squared()[0] for item in dict_to_merge[x]])]
      temp = [item.get_m1_squared()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'C':
    for x in xList:
      yList += [numpy.mean([item.get_C()[0] for item in dict_to_merge[x]])]
      temp = [item.get_C()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'sign':
    for x in xList:
      yList += [numpy.mean([item.get_sign()[0] for item in dict_to_merge[x]])]
      temp = [item.get_sign()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'sign_up':
    for x in xList:
      yList += [numpy.mean([item.get_sign_up()[0] for item in dict_to_merge[x]])]
      temp = [item.get_sign_up()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  elif variable == 'sign_down':
    for x in xList:
      yList += [numpy.mean([item.get_sign_down()[0] for item in dict_to_merge[x]])]
      temp = [item.get_sign_down()[1] ** 2 for item in dict_to_merge[x]]
      yErr += [math.sqrt(sum(temp)) / len(temp)]

  return xList, yList, yErr

