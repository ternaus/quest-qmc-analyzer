'''
example settings file
'''
from __future__ import division
import os

__author__ = 'Vladimir Iglovikov'

#folder in which different folders with different models output stored:
#ex:
# /home/valdimir/work/results
# In this folder I keep folders with different models outputs.
#ex: Lieb, kagome, square, etc
folder_with_different_models = os.path.join(os.path.expanduser("~"), 'work', 'results')