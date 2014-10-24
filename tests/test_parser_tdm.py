from __future__ import division
from unittest import TestCase
import os

import common.parser


__author__ = 'Vladimir Iglovikov'


class TestParser(TestCase):
  def setUp(self):
    dimension = 2
    self.time_indep_text = open(os.path.join(os.getcwd(), 'data', 'test.out')).read()
    self.time_dep_text = open(os.path.join(os.getcwd(), 'data', 'test.tdm.out')).read()
    self.geometry = open(os.path.join(os.getcwd(), 'data', 'test.geometry')).read()
    self.tparser = common.parser.Parser(self.time_indep_text, tdm=self.time_dep_text, geometry=self.geometry,
                                        dimension=dimension)

  def test_real_phase(self):
    # print self.tparser.get_ld_xx_real()
    for x, xe in self.tparser.get_ld_L().values():
      self.assertEquals(xe, 0)
    for x, xe in self.tparser.get_ld_T().values():
      self.assertEquals(xe, 0)
