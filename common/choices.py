from __future__ import division
import sys

__author__ = 'vladimir'

def list_choice(tlist_args):
  '''
  allows user to choose from the list, cases that he is interested in
  '''
  list_args = sorted(tlist_args)
  if len(list_args) == 0:
    print 'Zero variants. No choice'
    sys.exit(0)
  elif len(list_args) == 1:
      return [ list_args[0] ]
  else:
      for i, case in enumerate(list_args):
          print i,  ' - ',  case
      caseNumber = raw_input('Enter number of the case we are interested in\n')
      try:
          caseNumber = [int(caseNumber)]
      except:
          caseNumber = [int(j) for j in caseNumber.strip().split()]

      return [list_args[j] for j in caseNumber]
