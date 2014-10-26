from __future__ import division
import re

__author__ = 'Vladimir Iglovikov'


def extract_geometry(victim, **kwargs):
  parameter = kwargs['parameter']
  result = {}
  if parameter == 'coordinates':
    tm = re.search(r'(?<=Z)[0-9\s s\.]*', victim).group()
    for line in tm.strip().split("\n"):
      tp = line.strip().split()
      index, label, type, x, y, z = int(tp[0]), tp[1], int(tp[2]), float(tp[3]), float(tp[4]), float(tp[5])
      result[index] = (label, type, x, y, z)
    return result


def extract_tdm_data(victim, **kwargs):
  parameter = kwargs['parameter']

  if parameter == 'ld_xx_real':
    result = {}
    tm = re.findall(r'(?<=xx Current)([\s 0-9\.E+-]*)', victim.replace('+-', ''))
    for element in tm:
      key = element.strip().split('\n')[0].split()
      i = int(key[0])
      j = int(key[1])
      for line in element.strip().split('\n')[1:]:
        tl = line.strip().split()
        result[(i, j, float(tl[0]))] = (float(tl[1]), float(tl[2]))
    return result


def extract_non_tdm_data(victim, **kwargs):
  parameter = kwargs['parameter']

  if parameter == 'k_points':
    result = []
    tm = re.search(r'(?<=Class)[\s 0-9\.-]*', victim).group()
    for line in tm.strip().split('\n'):
      tl = line.strip().split()
      if len(tl) == 3:
        result += [(float(tl[1]), float(tl[2]))]
      elif len(tl) == 2:
        result += [(float(tl[0]), float(tl[1]))]
    return result
  elif parameter == "Mean Equal time Green's function":
    result = []
    tm = re.search(r"(?<=Mean Equal time Green's function:)[\s 0-9\.E+-]*", victim).group()
    tm = tm.replace('+-', '').strip().split("\n")
    for line in tm:
      tl = line.strip().split()
      orbit1, orbit2, x, y, z, symmetry, value, err = int(tl[0]), int(tl[1]), float(tl[2]), float(tl[3]), float(
        tl[4]), int(tl[5]), float(tl[6]), float(tl[7])
      result += [(orbit1, orbit2, x, y, z, symmetry, value, err)]
    return result

