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
    self.tparser = common.parser.Parser(self.time_indep_text, tdm=self.time_dep_text, dimension=dimension)

  def test_ld_xxp(self):
    print self.tparser.get_ld_xx_real()


