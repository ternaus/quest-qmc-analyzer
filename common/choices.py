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
    print 0, ' - all cases'
    for i, case in enumerate(list_args):
      print i + 1, ' - ', case
    caseNumber = raw_input('Enter number of the case we are interested in\n')
    try:
      caseNumber = [int(caseNumber)]
    except:
      caseNumber = [int(j) for j in caseNumber.strip().split()]
    else:
      if caseNumber[-1] == 0:
        return list_args

    result = [list_args[j - 1] for j in caseNumber]
    print result
    return [list_args[j - 1] for j in caseNumber]
