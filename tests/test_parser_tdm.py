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


  def test_L_match_T_at_mPi(self):
    self.assertTrue(abs(self.tparser.get_ld_L()[-3.14159][0] - self.tparser.get_ld_T()[-3.14159][0])) < (
    self.tparser.get_ld_L()[-3.14159][1] - self.tparser.get_ld_T()[-3.14159][1])


