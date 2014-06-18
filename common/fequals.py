from __future__ import division

__author__ = 'vladimir'

def equals(a, b, **kwargs):
  '''
  Function compares 2 float numbers.
  default precision is 5 decimal places.
  @param a:
  @param b:
  @param kwargs:
    place - number of places after the decimal point that we
    want these numbers being equal to
  @return: True if |a - b| < epsilon, False otherwise
  '''
  if 'places' in kwargs:
    epsilon = 1.0 / 10**kwargs['places']
  else:
    epsilon = 1e-5
  return abs(a-b) < epsilon