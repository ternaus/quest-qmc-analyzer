from __future__ import division

from unittest import TestCase
import os

__author__ = 'Vladimir Iglovikov'


class TestCommon(TestCase):
  def test_extract_tdm_SxSx(self):
    time_indep_text = open(os.path.join(os.getcwd(), 'data', 'Lieb_1385474159.53.out')).read()
    time_dep_text = open(os.path.join(os.getcwd(), 'data', 'Lieb_1385474159.53.tdm.out')).read()

    parameter = 'XX Spin correlation function:'


