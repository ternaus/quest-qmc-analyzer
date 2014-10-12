from __future__ import division

__author__ = 'Vladimir Iglovikov'
import statsmodels.api as sm


def approx(x_list, y_list, y_err):
  X = sm.add_constant(x_list)
  wls_fit = sm.WLS(y_list, X, weights=[1.0 / tx for tx in y_err]).fit()
  # wls_fit = sm.WLS(y_list, X).fit()

  c = wls_fit.params[0]
  m = wls_fit.params[1]
  std_err = wls_fit.bse[0]

  approxx = [0] + x_list

  approxy = [m * tt + c for tt in approxx]

  return approxx, approxy, m, c, std_err