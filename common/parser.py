#!/usr/bin/env python
'''
This class get's data from the file in which QUEST saves all data
and saves it as a fileds of the class

'''
from __future__ import division
import re
import math
import sys

import common
import common.extract_data


__author__ = 'vladimir'


class Parser:
  def __init__(self, fileText, **kwargs):
    self.fileText = fileText
    self.dimension = kwargs['dimension']

    if 'tdm' not in kwargs or kwargs['tdm'] == '':
      self.fileText_tdm = ''
    else:
      self.fileText_tdm = open(kwargs['tdm']).read()

    if 'geometry' not in kwargs or kwargs['geometry'] == '':
      self.geometry = ''
    else:
      self.geometry = open(kwargs['geometry']).read()

    self.u = None
    self.num_slices = None  # Number of time slices
    self.t_up = None
    self.t_down = None
    self.mu_up = None
    self.mu_down = None
    self.nSites = None
    self.beta = None
    self.energy = None
    self.n_up = None
    self.n_down = None
    self.rho = None
    self.struct_XX_F = None
    self.struct_XX_AF = None
    self.SxSx = None
    self.SzSz = None
    self.SxSx_tdm = None
    self.SzSz_tdm = None
    self.X_F = None
    self.dtau = None
    self.DO = None  # Double occupancy
    self.s_wave = None  # s-wave pairing
    self.energy_hop = None  # hopping term
    self.L = None  # Effective number of unit cells in one direction
    self.nx = None  # Number of unit cells in x direction
    self.ny = None  # Number of unit cells in y direction
    self.pairing = None  # Pairing correlation function.
    self.pairing_vs_x = None  # Pairing vs x
    self.orbits = None  # sorted list of the orbits
    self.num_orbits = None  # number of orbits
    self.global_sites = None  # Number of global move sites
    self.global_accept = None  # Global move accept rate
    self.m2 = None  # Square of the magnetisation
    self.m = None  # Magnetisation
    self.bc = None  # Binder cumulant, defined as self.m2 / (self.m)^2
    self.n_up2 = None  # Square of the upspin density
    self.n_down2 = None  # Square of the downspin density
    self.density_correlation_upup = None  # Density up -density up correlation function
    self.density_correlation_updn = None  # Density up -density down correlation function
    self.density_correlation = None  # Density - density correlation function
    self.n0 = None  # Density on the 0th orbital
    self.n1 = None  # Density on the 1th orbital
    self.green = None  # Mean Equal time green function
    self.green_up = None  # Up equal time green function
    self.green_down = None  # Down equal time green function
    self.n01_n0_n1 = None  # <n0 n1> - <n0> <n1>
    self.n12_n1_n2 = None  # <n1 n2> - <n1> <n2>
    self.m0_squared = None  # Square of the magnetisation on the orbital 0
    self.m1_squared = None  # Square of the magnetisation on the orbital 1
    self.C = None  # Specific heat
    self.FT_pairing = None  # Fourier transfrom of the pairing correlation function
    self.FT_pairing_00 = None  # Fourier transform of the pairing correlation function for the 00 orbitals
    self.FT_pairing_10 = None  # Fourier transform of the pairing correlation function for the 10 orbitals
    self.FT_pairing_11 = None  # Fourier transform of the pairing correlation function for the 11 orbitals
    self.FT_pairing_21 = None  # Fourier transform of the pairing correlation function for the 21 orbitals
    self.kx_points = None  # k points in the x direction
    self.ky_points = None  # k points in the y direction
    self.sign = None  # average sign of the determinant
    self.sign_up = None  # average sign for the determinant corresponding to the spin up electrons
    self.sign_down = None  # average sign for the determinant corresponding to the spin down electrons
    self.sign_normalized = None  # <S> - <S_up> <S_dn>
    self.sign_up_down = None  # sign_up * sign_down
    self.ld_xx_real = None  # real space current current correlation function
    self.k_points = None  # list of tuples of k points
    self.ld_T = None  # transverse response
    self.ld_L = None  # longitudinal response
    self.coordinates = None  # coordinates of the sites
    self.tau_list = None  # list of tau
    self.kx = None  # Hopping energy measured along x axis
    self.ld_w = None  #ld_xx at q = 0, w != 0

  def get_num_slices(self):
    if self.num_slices == None:
      self.num_slices = int(re.search('(?<=Time slice - L :)\s+.?\d+', self.fileText).group(0))
    return self.num_slices

  def get_ld_w(self):
    # TODO after site numeration is fixed to remove the hack.
    if self.ld_w == None:
      self.ld_w = {}
      for omega in range(self.get_num_slices()):
        tp_real = 0
        tp_im = 0
        tp_real_err = 0
        for site1, site2, tau in self.get_ld_xx_real():
          tp_real += math.cos(2 * math.pi * omega * tau / self.get_beta()) * self.get_ld_xx_real()[(site1, site2, tau)][
            0]
          tp_real_err += math.cos(2 * math.pi * omega * tau / self.get_beta()) * \
                         self.get_ld_xx_real()[(site1, site2, tau)][1]
          tp_im += math.sin(2 * math.pi * omega * tau / self.get_beta()) * self.get_ld_xx_real()[(site1, site2, tau)][0]

        if abs(tp_im) > 0.02:
          print 'ld_w error'
          print 'omega = ', omega
          print 'tp_im = ', tp_im
          sys.exit(0)

        self.ld_w[omega] = (
          self.get_dtau() * tp_real / self.get_nSites(), self.get_dtau() * tp_real_err / self.get_nSites())
    return self.ld_w

  def get_kx(self):
    if self.kx == None:
      # find minimum nonzero x separation
      separation_list = list(set([tx[2] for tx in self.get_green()]))
      separation_list.sort()
      min_x = separation_list[1]
      for value in self.get_green():
        if value[2] == min_x and value[3] == 0 and value[4] == 0:
          self.kx = (value[6], value[7])
    return self.kx


  def get_coordinates(self):
    if self.coordinates == None:
      self.coordinates = common.extract_data.extract_geometry(self.geometry, parameter='coordinates')
    return self.coordinates

  def get_tau_list(self):
    if self.tau_list == None:
      self.tau_list = list(set([tx[2] for tx in self.get_ld_xx_real().keys()]))
      self.tau_list.sort()
    return self.tau_list

  def get_k_points(self):
    if self.k_points == None:
      self.k_points = common.extract_data.extract_non_tdm_data(self.fileText, parameter='k_points',
                                                               dimension=self.dimension)

      self.kx_points = list(set([tx[0] for tx in self.get_k_points()]))

      self.kx_points.sort()
      self.ky_points = list(set([tx[1] for tx in self.get_k_points()]))

      self.ky_points.sort()
    return self.k_points


  def get_ld_L(self):
    # TODO after site numeration is fixed to remove the hack.
    if self.ld_L == None:
      self.ld_L = {}
      for qx in self.get_kx_points():
        if qx > 0:
          continue
        tp_real = 0
        tp_im = 0
        tp_real_err = 0
        for site1, site2, tau in self.get_ld_xx_real():
          lx1, ly1 = self.get_coordinates()[site1 - 1][2], self.get_coordinates()[site1 - 1][3]
          lx2, ly2 = self.get_coordinates()[site2 - 1][2], self.get_coordinates()[site2 - 1][3]

          m_cos = math.cos((lx1 - lx2) * qx)
          tp_real += m_cos * self.get_ld_xx_real()[(site1, site2, tau)][0]
          tp_real_err += m_cos * self.get_ld_xx_real()[(site1, site2, tau)][1]
          tp_im += math.sin((lx1 - lx2) * qx) * self.get_ld_xx_real()[(site1, site2, tau)][0]

        if abs(tp_im) > 0.02:
          print 'ld_L error'
          print 'qx = ', qx
          print 'tp_im = ', tp_im
          sys.exit(0)

        self.ld_L[qx] = (
          self.get_dtau() * tp_real / self.get_nSites(), self.get_dtau() * tp_real_err / self.get_nSites())
    return self.ld_L

  def get_ld_T(self):
    '''

    @return: dictionary where key is qy and value is real and imaginary part
    '''
    # TODO after site numeration is fixed to remove the hack.
    if self.ld_T == None:
      self.ld_T = {}

      for qy in self.get_ky_points()[:4]:
        # if qy >= 0 or qy == -3.14159:
        # continue
        tp_real = 0
        tp_im = 0
        tp_real_err = 0
        for site1, site2, tau in self.get_ld_xx_real():
          lx1, ly1 = self.get_coordinates()[site1 - 1][2], self.get_coordinates()[site1 - 1][3]
          lx2, ly2 = self.get_coordinates()[site2 - 1][2], self.get_coordinates()[site2 - 1][3]
          m_cos = math.cos((ly1 - ly2) * qy)
          tp_real += m_cos * self.get_ld_xx_real()[(site1, site2, tau)][0]
          tp_real_err += m_cos * self.get_ld_xx_real()[(site1, site2, tau)][1]

          tp_im += math.sin((ly1 - ly2) * qy) * self.get_ld_xx_real()[(site1, site2, tau)][0]

        if abs(tp_im) > 1e-5:
          print 'ld_T error'
          print 'qy = ', qy
          print 'tp_im = ', tp_im
          sys.exit(0)

        self.ld_T[qy] = (
          self.get_dtau() * tp_real / self.get_nSites(), self.get_dtau() * tp_real_err / self.get_nSites())
    return self.ld_T


  def get_ld_xx_real(self):
    '''
    @return: dictionary, where key index of site 1, index of site 2, tau.
     value - value of current xx + errorbars.
    '''
    if self.ld_xx_real == None:
      self.ld_xx_real = common.extract_data.extract_tdm_data(self.fileText_tdm, parameter='ld_xx_real')
    return self.ld_xx_real

  def get_num_orbits(self):
    if self.num_orbits == None:
      self.num_orbits = len(self.get_orbits())
    return self.num_orbits

  def get_u(self):
    if self.u == None:
      self.u = float(re.search('(?<=U :)\s+.?\d+\.\d+', self.fileText).group(0))
    return self.u

  def get_t_up(self):
    if self.t_up == None:
      self.t_up = re.search('t_up :.*', self.fileText).group(0).replace('t_up :', '').strip().split()
      self.t_up = [float(tx) for tx in self.t_up]

    return self.t_up

  def get_t_down(self):
    if self.t_down == None:
      self.t_down = re.search('t_dn :.*', self.fileText).group(0).replace('t_dn :', '').strip().split()
      self.t_down = [float(tx) for tx in self.t_down]
    return self.t_down

  def get_mu_up(self):
    if self.mu_up == None:
      self.mu_up = float(re.search('(?<=mu_up :)\s+.?\d+\.\d+', self.fileText).group(0))
    return self.mu_up

  def get_mu_down(self):
    if self.mu_down == None:
      self.mu_down = float(re.search('(?<=mu_dn :)\s+.?\d+\.\d+', self.fileText).group(0))
    return self.mu_down

  def get_nSites(self):
    if self.nSites == None:
      self.nSites = int(re.search('(?<=Number of sites :)\s+.?\d+', self.fileText).group(0))

    return self.nSites

  def get_global_sites(self):
    if self.global_sites == None:
      self.global_sites = int(re.search('(?<=Global move number of sites :)\s+.?\d+', self.fileText).group(0))
    return self.global_sites

  def get_global_accept(self):
    if self.get_global_sites() > 0:
      if self.global_accept == None:
        self.global_accept = float(re.search('(?<=Global move accept rate :)\s+.?\d+', self.fileText).group(0))
    return self.global_accept

  def get_L(self):
    if self.L == None:
      self.L = math.sqrt(self.get_nx() * self.get_ny())
    return self.L

  def get_nx(self):
    if self.nx == None:
      self.nx = len(self.get_kx_points())
    return self.nx

  def get_ny(self):
    if self.ny == None:
      self.ny = len(self.get_ky_points())
    return self.ny

  def get_beta(self):
    if self.beta == None:
      self.beta = float(re.search('(?<=beta :)\s+.?\d+\.\d+', self.fileText).group(0))
    return self.beta

  def get_n_up(self):
    if self.n_up == None:
      m = re.search('(?<=Up spin occupancy :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.n_up = (float(m.groups()[0]), float(m.groups()[1]))
    return self.n_up

  def get_n_up2(self):
    if self.n_up2 == None:
      self.n_up2 = (self.get_n_up()[0] ** 2, 2 * math.sqrt(2) * self.get_n_up()[1] * self.get_n_up()[0])
    return self.n_up2

  def get_n_down(self):
    if self.n_down == None:
      m = re.search('(?<=Down spin occupancy :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)',
                    self.fileText)
      self.n_down = (float(m.groups()[0]), float(m.groups()[1]))
    return self.n_down

  def get_n_down2(self):
    if self.n_down2 == None:
      self.n_down2 = (self.get_n_down()[0] ** 2, 2 * math.sqrt(2) * self.get_n_down()[1] * self.get_n_down()[0])
    return self.n_down2

  def get_DO(self):
    if self.DO == None:
      try:
        m = re.search('(?<=Double occupancy :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
        self.DO = (float(m.groups()[0]), float(m.groups()[1]))
      except:
        print 'here'
        if self.get_u() == 0:
          self.DO = (self.get_n_up()[0] * self.get_n_down()[0], math.sqrt(
            (self.get_n_down()[1] * self.get_n_up()[0]) ** 2 + (self.get_n_down()[0] * self.get_n_up()[1]) ** 2))
        else:
          m = re.search('(?<= <U\*N_up\*N_dn\>  :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)',
                        self.fileText)

          self.DO = (float(m.groups()[0]) / self.get_u(), float(m.groups()[1]) / self.get_u())
    return self.DO

  def get_m2(self):
    if self.m2 == None:
      try:
        m = re.search('(?<=Magnetizatiion squared :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)',
                      self.fileText)
        self.m2 = (float(m.groups()[0]), float(m.groups()[1]))
      except:
        self.m2 = (self.get_n_up()[0] + self.get_n_down()[0] - 2 * self.get_DO()[0],
                   math.sqrt(self.get_n_up()[1] ** 2 + self.get_n_down()[1] ** 2 + 4 * self.get_DO()[1] ** 2))
    return self.m2

  def get_m(self):
    if self.m == None:
      self.m = (
        self.get_n_up()[0] - self.get_n_down()[0], math.sqrt(self.get_n_up()[1] ** 2 + self.get_n_down()[1] ** 2))
    return self.m

  def get_bc(self):
    if self.bc == None:
      self.bc = (self.get_m2()[0] / self.get_m()[0] ** 2, math.sqrt((self.get_m2()[1] / self.get_m()[0]) ** 2 + (
        2 * self.get_m2()[0] * self.get_m()[1] / self.get_m()[0] ** 3) ** 2))
    return self.bc

  def get_rho(self):
    if self.rho == None:
      m = re.search('(?<=Density :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.rho = (float(m.groups()[0]), float(m.groups()[1]))
      if round(self.get_n_up()[0] + self.get_n_down()[0] - self.rho[0], 6) != 0:
        print self.get_n_up(), self.get_n_down(), self.rho
      assert round(self.get_n_up()[0] + self.get_n_down()[0] - self.rho[0], 6) == 0
    return self.rho

  def get_C(self):
    if self.C == None:
      m = re.search('(?<=Specific heat :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.C = (float(m.groups()[0]), float(m.groups()[1]))
    return self.C

  def get_energy_hop(self):
    if self.energy_hop == None:
      m = re.search('(?<=Hopping energy :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.energy_hop = (float(m.groups()[0]), float(m.groups()[1]))
    return self.energy_hop

  def get_energy(self):
    if self.energy == None:
      try:
        m = re.search('(?<=Total [e, E]nergy :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)',
                      self.fileText)
        self.energy = (float(m.groups()[0]), float(m.groups()[1]))
      except:
        print self.fileText
    return self.energy

  def get_sign(self):
    if self.sign == None:
      m = re.search('(?<=Avg sign :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.sign = (float(m.groups()[0]), float(m.groups()[1]))
    return self.sign

  def get_sign_up(self):
    if self.sign_up == None:
      m = re.search('(?<=Avg up sign :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.sign_up = (float(m.groups()[0]), float(m.groups()[1]))
    return self.sign_up

  def get_sign_down(self):
    if self.sign_down == None:
      m = re.search('(?<=Avg dn sign :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.sign_down = (float(m.groups()[0]), float(m.groups()[1]))
    return self.sign_down

  def get_sign_normalized(self):
    if self.sign_normalized == None:
      self.sign_normalized = (self.get_sign()[0] - self.get_sign_up()[0] * self.get_sign_down()[0],
                              math.sqrt(
                                self.get_sign()[1] ** 2 + self.get_sign_down()[0] ** 2 * self.get_sign_up()[1] ** 2 +
                                self.get_sign_up()[0] ** 2 * self.get_sign_down()[1] ** 2))
    return self.sign_normalized

  def get_sign_up_down(self):
    if self.sign_up_down == None:
      self.sign_up_down = (self.get_sign_up()[0] * self.get_sign_down()[0],
                           math.sqrt(self.get_sign_down()[0] ** 2 * self.get_sign_up()[1] ** 2 +
                                     self.get_sign_up()[0] ** 2 * self.get_sign_down()[1] ** 2))
    return self.sign_up_down


  def get_struct_XX_F(self):
    if self.struct_XX_F == None:
      m = re.search('(?<=XX Ferro structure factor :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)',
                    self.fileText)
      self.struct_XX_F = (float(m.groups()[0]), float(m.groups()[1]))
    return self.struct_XX_F

  def get_struct_XX_AF(self):
    if self.struct_XX_AF == None:
      m = re.search('(?<=XX AF structure factor :)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)',
                    self.fileText)
      self.struct_XX_AF = (float(m.groups()[0]), float(m.groups()[1]))
    return self.struct_XX_AF

  def get_SxSx(self):
    if self.SxSx == None:
      self.SxSx = common.extract_data.extract_non_tdm_data(self.fileText, parameter='XX Spin correlation function:')
    return self.SxSx

  def get_SzSz(self):
    if self.SzSz == None:
      self.SzSz = common.extract_data.extract_non_tdm_data(self.fileText, parameter='ZZ Spin correlation function:')
    return self.SzSz

  def get_pairing(self):
    if self.pairing == None:
      self.pairing = common.extract_data.extract_non_tdm_data(self.fileText, parameter='Pairing correlation function:',
                                                              dimension=self.dimension)
    return self.pairing

  def get_SxSx_tdm(self):
    if self.SxSx_tdm == None:
      self.SxSx_tdm = common.extract_data.extract_tdm_data(self.fileText_tdm, parameter='SxSx')
    return self.SxSx_tdm

  def get_SzSz_tdm(self):
    if self.SzSz_tdm == None:
      self.SzSz_tdm = common.extract_data.extract_tdm_data(self.fileText_tdm, parameter='SzSz')
    return self.SzSz_tdm

  def get_dtau(self):
    if self.dtau == None:
      self.dtau = float(re.search('(?<=dtau :)\s+.?\d+\.\d+', self.fileText).group(0))
    return self.dtau


  def get_s_wave(self):
    if self.s_wave == None:
      m = re.search('(?<=s-wave)\s+(\-?\d+\.\d+E?-?\+?\d+)\s+\+?-?\s+(\d+\.\d+E?\+?-?\d+)', self.fileText)
      self.s_wave = (float(m.groups()[0]), float(m.groups()[1]))
    return self.s_wave


  def get_kx_points(self):
    if self.kx_points == None:
      self.get_k_points()
    return self.kx_points

  def get_ky_points(self):
    if self.ky_points == None:
      self.get_k_points()
    return self.ky_points

  def get_pairing_vs_x(self):
    if self.pairing_vs_x == None:
      resultDict = {}
      for key, value in self.get_pairing().iteritems():
        if key[3] == 0:  # we are interested only in values where y jump is 0
          newKey = (key[0], key[1])
          if newKey not in resultDict:
            resultDict[newKey] = [(key[2], value[0], value[1])]
          else:
            resultDict[newKey] += [(key[2], value[0], value[1])]

      for key in resultDict:
        temp = resultDict[key]
        temp.sort()
        resultDict[key] = [[tx[0] for tx in temp], [tx[1] for tx in temp], [tx[2] for tx in temp]]
      self.pairing_vs_x = resultDict
    return self.pairing_vs_x

  def get_orbits(self):
    if self.orbits == None:
      self.orbits = set()
      for key in self.get_pairing():
        self.orbits.add(key[0])
        self.orbits.add(key[1])
      self.orbits = list(self.orbits)
      self.orbits.sort()
    return self.orbits

  def get_density_correlation_up_up(self):
    if self.density_correlation_upup == None:
      self.density_correlation_upup = common.extract_data.extract_non_tdm_data(self.fileText,
                                                                               parameter='Density-density correlation fn: (up-up)')
    return self.density_correlation_upup

  def get_density_correlation_up_down(self):
    if self.density_correlation_updn == None:
      self.density_correlation_updn = common.extract_data.extract_non_tdm_data(self.fileText,
                                                                               parameter='Density-density correlation fn: (up-dn)')
    return self.density_correlation_updn

  def get_density_correlation(self):
    if self.density_correlation == None:
      self.density_correlation = {}
      for key in self.get_density_correlation_up_down().keys():
        value_up_up = self.get_density_correlation_up_up()[key]
        value_up_down = self.get_density_correlation_up_down()[key]
        self.density_correlation[key] = (2 * (value_up_up[0] + value_up_down[0]),
                                         2 * math.sqrt(value_up_up[1] ** 2 + value_up_down[1] ** 2))
    return self.density_correlation

  def get_green(self):
    if self.green == None:
      self.green = common.extract_data.extract_non_tdm_data(self.fileText,
                                                            parameter="Mean Equal time Green's function")
    return self.green

  def get_green_up(self):
    if self.green_up == None:
      self.green_up = common.extract_data.extract_non_tdm_data(self.fileText,
                                                               parameter="Up Equal time Green's function:")
    return self.green_up

  def get_green_down(self):
    if self.green_down == None:
      self.green_down = common.extract_data.extract_non_tdm_data(self.fileText,
                                                                 parameter="Down Equal time Green's function:")
    return self.green_down

  def get_n0(self):
    if self.n0 == None:
      self.n0 = ((2 - self.get_green_up()[(0, 0, 0, 0, 0)][0] - self.get_green_down()[(0, 0, 0, 0, 0)][0]),
                 math.sqrt(
                   self.get_green_up()[(0, 0, 0, 0, 0)][1] ** 2 + self.get_green_down()[(0, 0, 0, 0, 0)][1] ** 2))
    return self.n0

  def get_n1(self):
    if self.n1 == None:
      self.n1 = ((2 - self.get_green_up()[(1, 1, 0, 0, 0)][0] - self.get_green_down()[(1, 1, 0, 0, 0)][0]),
                 math.sqrt(
                   self.get_green_up()[(1, 1, 0, 0, 0)][1] ** 2 + self.get_green_down()[(1, 1, 0, 0, 0)][1] ** 2))
    return self.n1

  def get_n01_normalized(self):
    if self.n01_n0_n1 == None:
      self.n01_n0_n1 = (self.get_density_correlation()[0, 1, 0.5, 0, 0][0] - self.get_n0()[0] * self.get_n1()[0],
                        math.sqrt(self.get_density_correlation()[0, 1, 0.5, 0, 0][1] ** 2 + (
                          self.get_n0()[0] * self.get_n1()[1]) ** 2 + (self.get_n1()[0] * self.get_n0()[1]) ** 2))
    return self.n01_n0_n1

  def get_n12_normalized(self):
    if self.n12_n1_n2 == None:
      self.n12_n1_n2 = (self.get_density_correlation()[1, 2, -0.5, 0.5, 0][0] - self.get_n1()[0] ** 2, math.sqrt(
        self.get_density_correlation()[1, 2, -0.5, 0.5, 0][1] ** 2 + 2 * (self.get_n1()[0] * self.get_n1()[1]) ** 2))
    return self.n12_n1_n2

  def get_m0_squared(self):
    if self.m0_squared == None:
      # TODO errobars
      self.m0_squared = (self.get_n0()[0] - 2 * self.get_density_correlation_up_down()[(0, 0, 0, 0, 0)][0], 0)
    return self.m0_squared

  def get_m1_squared(self):
    if self.m1_squared == None:
      # TODO errobars
      self.m1_squared = (self.get_n1()[0] - 2 * self.get_density_correlation_up_down()[(1, 1, 0, 0, 0)][0], 0)
    return self.m1_squared