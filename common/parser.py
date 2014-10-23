#!/usr/bin/env python
'''
This class get's data from the file in which QUEST saves all data
and saves it as a fileds of the class

'''
from __future__ import division
import re
import math
import numpy

import common
import common.extract_data


__author__ = 'vladimir'


class Parser:
  def __init__(self, fileText, **kwargs):
    self.fileText = fileText
    self.dimension = kwargs['dimension']

    if 'tdm' in kwargs:
      self.fileText_tdm = kwargs['tdm']

    self.u = None
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
    self.DO = None #Double occupancy
    self.s_wave = None #s-wave pairing
    self.energy_hop = None #hopping term
    self.L = None #Effective number of unit cells in one direction
    self.nx = None #Number of unit cells in x direction
    self.ny = None #Number of unit cells in y direction
    self.pairing = None #Pairing correlation function.
    self.pairing_vs_x = None #Pairing vs x
    self.orbits = None #sorted list of the orbits
    self.num_orbits = None #number of orbits
    self.global_sites = None #Number of global move sites
    self.global_accept = None #Global move accept rate
    self.m2 = None #Square of the magnetisation
    self.m = None #Magnetisation
    self.bc = None #Binder cumulant, defined as self.m2 / (self.m)^2
    self.n_up2 = None #Square of the upspin density
    self.n_down2 = None #Square of the downspin density
    self.density_correlation_upup = None # Density up -density up correlation function
    self.density_correlation_updn = None # Density up -density down correlation function
    self.density_correlation = None # Density - density correlation function
    self.n0 = None #Density on the 0th orbital
    self.n1 = None #Density on the 1th orbital
    self.green = None #Mean Equal time green function
    self.green_up = None # Up equal time green function
    self.green_down = None # Down equal time green function
    self.n01_n0_n1 = None #<n0 n1> - <n0> <n1>
    self.n12_n1_n2 = None #<n1 n2> - <n1> <n2>
    self.m0_squared = None #Square of the magnetisation on the orbital 0
    self.m1_squared = None #Square of the magnetisation on the orbital 1
    self.C = None #Specific heat
    self.k_grid = None #Grid of the k points corresponding to the different classes
    self.FT_pairing = None #Fourier transfrom of the pairing correlation function
    self.FT_pairing_00 = None #Fourier transform of the pairing correlation function for the 00 orbitals
    self.FT_pairing_10 = None #Fourier transform of the pairing correlation function for the 10 orbitals
    self.FT_pairing_11 = None #Fourier transform of the pairing correlation function for the 11 orbitals
    self.FT_pairing_21 = None #Fourier transform of the pairing correlation function for the 21 orbitals
    self.kx_points = None #k points in the x direction
    self.ky_points = None #k points in the y direction
    self.sign = None  # average sign of the determinant
    self.sign_up = None  # average sign for the determinant corresponding to the spin up electrons
    self.sign_down = None  #average sign for the determinant corresponding to the spin down electrons
    self.sign_normalized = None  # <S> - <S_up> <S_dn>
    self.sign_up_down = None  # sign_up * sign_down
    self.ld_xx_real = None  # real space current current correlation function
    self.k_points = None  # list of tuples of k points
    self.ld_T = None  # transverse response
    self.ld_L = None  # longitudinal response
    self.coordinates = None  # coordinates of the sites
    self.tau_list = None  # list of tau
    self.rho_s = None  # rho_s

  def get_tau_list(self):
    if self.tau_list == None:
      self.tau_list = list(set([tx[2] for tx in self.get_ld_xx_real().keys()]))
      self.tau_list.sort()
    return self.tau_list

  def get_k_points(self):
    if self.k_points == None:
      self.k_points = common.extract_data.extract_non_tdm_data(self.fileText, parameter='k_points',
                                                               dimension=self.dimension)
    return self.k_points

  def get_ld_T(self):
    if self.ld_T == None:
      result = {}
      for lx, ly in self.get_coordinates():
        for tau in self.get_tau_list():
          for qy in self.get_ky_points():
            result[qy] = math.cos(ly * qy) * self.get_ld_xx_real()[(lx, ly, tau)][0]

    return self.ld_T

  def get_ld_L(self):
    if self.ld_L == None:
      result = {}
      for lx, ly in self.get_coordinates():
        for tau in self.get_tau_list():
          for qx in self.get_kx_points():
            result[qx] = math.cos(lx * qx) * self.get_ld_xx_real()[(lx, ly, tau)][0]

    return self.ld_L

  # TODO join get_ld_L and get_ld_T.

  def get_rho_s(self):
    if self.rho_s == None:
      self.rho_s = 0.25 * (self.get_ld_L() - self.get_ld_T())
    return self.rho_s

  def get_ld_xx_real(self):
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
      self.t_down = re.search('t_dn :.*', self.fileText).group(0).replace('t_up :', '').strip().split()
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
      if self.dimension == 1:
        assert (self.nSites == self.get_nx())

        # if self.dimension == 2:
        # try:
        #       assert (self.nSites == self.get_num_orbits() * self.get_nx() * self.get_ny())
        #   except:
        #       self.nx = self.ny = int(math.sqrt(self.nSites / 2.0))

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
      self.pairing = common.extract_data.extract_non_tdm_data(self.fileText, parameter='Pairing correlation function:', dimension=self.dimension)
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

  def get_FT_pairing(self):
    if self.FT_pairing == None:
      self.FT_pairing = common.extract_data.extract_non_tdm_data(self.fileText,
                                                                 parameter='FT of Pairing correlation fn:', FT=True)
    return self.FT_pairing

  def get_FT_pairing_00(self):
    if self.FT_pairing_00 == None:
      self.FT_pairing_00 = numpy.zeros((self.get_nx(), self.get_ny()))
      for kx in range(self.get_nx()):
        for ky in range(self.get_ny()):
          self.FT_pairing_00[(kx, ky)] = -1

      temp = 0
      for value in self.get_FT_pairing():
        if value[1] == 0 and value[2] == 0:
          for kx, ky in self.get_k_grid()[value[0]]:
            temp += 1
            if self.get_nx() % 2 == 0:
              kx_index = round((kx + math.pi) * self.get_nx() / (2 * math.pi))
            elif self.get_nx() % 2 == 1:
              kx_index = round(((kx + math.pi) * self.get_nx() / (math.pi) - 1) / 2)

            if self.get_ny() % 2 == 0:
              ky_index = round((ky + math.pi) * self.get_ny() / (2 * math.pi))
            elif self.get_ny() % 2 == 1:
              ky_index = round(((ky + math.pi) * self.get_ny() / (math.pi) - 1) / 2)

            self.FT_pairing_00[(kx_index, ky_index)] = value[3]

      assert temp == self.get_nx() * self.get_ny() #Check if we cover all points

    return self.FT_pairing_00

  def get_FT_pairing_10(self):
    if self.FT_pairing_10 == None:
      self.FT_pairing_10 = numpy.zeros((self.get_nx(), self.get_ny()))
      for kx in range(self.get_nx()):
        for ky in range(self.get_ny()):
          self.FT_pairing_10[(kx, ky)] = -1

      temp = 0
      for value in self.get_FT_pairing():
        if value[1] == 1 and value[2] == 0:
          for kx, ky in self.get_k_grid()[value[0]]:
            temp += 1
            if self.get_nx() % 2 == 0:
              kx_index = round((kx + math.pi) * self.get_nx() / (2 * math.pi))
            elif self.get_nx() % 2 == 1:
              kx_index = round(((kx + math.pi) * self.get_nx() / (math.pi) - 1) / 2)

            if self.get_ny() % 2 == 0:
              ky_index = round((ky + math.pi) * self.get_ny() / (2 * math.pi))
            elif self.get_ny() % 2 == 1:
              ky_index = round(((ky + math.pi) * self.get_ny() / (math.pi) - 1) / 2)

            self.FT_pairing_10[(kx_index, ky_index)] = value[3]

      assert temp == self.get_nx() * self.get_ny() #Check if we cover all points

    return self.FT_pairing_10

  def get_FT_pairing_11(self):
    if self.FT_pairing_11 == None:
      self.FT_pairing_11 = numpy.zeros((self.get_nx(), self.get_ny()))
      for kx in range(self.get_nx()):
        for ky in range(self.get_ny()):
          self.FT_pairing_11[(kx, ky)] = -1

      temp = 0
      for value in self.get_FT_pairing():
        if value[1] == 1 and value[2] == 1:
          for kx, ky in self.get_k_grid()[value[0]]:
            temp += 1
            if self.get_nx() % 2 == 0:
              kx_index = round((kx + math.pi) * self.get_nx() / (2 * math.pi))
            elif self.get_nx() % 2 == 1:
              kx_index = round(((kx + math.pi) * self.get_nx() / (math.pi) - 1) / 2)

            if self.get_ny() % 2 == 0:
              ky_index = round((ky + math.pi) * self.get_ny() / (2 * math.pi))
            elif self.get_ny() % 2 == 1:
              ky_index = round(((ky + math.pi) * self.get_ny() / (math.pi) - 1) / 2)

            self.FT_pairing_11[(kx_index, ky_index)] = value[3]

      assert temp == self.get_nx() * self.get_ny() #Check if we cover all points

    return self.FT_pairing_11

  def get_FT_pairing_21(self):
    if self.FT_pairing_21 == None:
      self.FT_pairing_21 = numpy.zeros((self.get_nx(), self.get_ny()))
      for kx in range(self.get_nx()):
        for ky in range(self.get_ny()):
          self.FT_pairing_21[(kx, ky)] = -1

      temp = 0
      for value in self.get_FT_pairing():
        if value[1] == 2 and value[2] == 1:
          for kx, ky in self.get_k_grid()[value[0]]:
            temp += 1
            if self.get_nx() % 2 == 0:
              kx_index = round((kx + math.pi) * self.get_nx() / (2 * math.pi))
            elif self.get_nx() % 2 == 1:
              kx_index = round(((kx + math.pi) * self.get_nx() / (math.pi) - 1) / 2)

            if self.get_ny() % 2 == 0:
              ky_index = round((ky + math.pi) * self.get_ny() / (2 * math.pi))
            elif self.get_ny() % 2 == 1:
              ky_index = round(((ky + math.pi) * self.get_ny() / (math.pi) - 1) / 2)

            self.FT_pairing_21[(kx_index, ky_index)] = value[3]

      assert temp == self.get_nx() * self.get_ny() #Check if we cover all points

    return self.FT_pairing_21

  def get_k_grid(self):
    if self.k_grid == None:
      self.k_grid = common.extract_data.extract_non_tdm_data(self.fileText,
                                                             parameter="Grid for Green's function", k_grid=True, dimension=self.dimension)
      grid_points_x = set()
      grid_points_y = set()
      temp = 0 #We need to check if number of k points in the k_grid is equal to the number of unit cells
      for item in self.k_grid.values():
        temp += len(item)
        if self.dimension == 2:
          for (kx, ky) in item:
            grid_points_x.add(kx)
            grid_points_y.add(ky)
        elif self.dimension == 1:
          for kx in item:
            grid_points_x.add(kx)
      self.kx_points = list(grid_points_x)
      self.ky_points = list(grid_points_y)
      self.kx_points.sort()
      self.ky_points.sort()

    return self.k_grid

  def get_kx_points(self):
    if self.kx_points == None:
      self.get_k_grid()
    return self.kx_points

  def get_ky_points(self):
    if self.ky_points == None:
      self.get_k_grid()
    return self.ky_points

  def get_pairing_vs_x(self):
    if self.pairing_vs_x == None:
      resultDict = {}
      for key, value in self.get_pairing().iteritems():
        if key[3] == 0:#we are interested only in values where y jump is 0
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
                                                            parameter="Mean Equal time Green's function:")
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
      #TODO errobars
      self.m0_squared = (self.get_n0()[0] - 2 * self.get_density_correlation_up_down()[(0, 0, 0, 0, 0)][0], 0)
    return self.m0_squared

  def get_m1_squared(self):
    if self.m1_squared == None:
      #TODO errobars
      self.m1_squared = (self.get_n1()[0] - 2 * self.get_density_correlation_up_down()[(1, 1, 0, 0, 0)][0], 0)
    return self.m1_squared