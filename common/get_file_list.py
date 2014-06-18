from __future__ import division
import os

import parser

__author__ = 'vladimir'


def get_filelist(modelName, path):
  """
  Returns list of the files with time independent measurements
  """
  #get list of the filenames that have information that we need
  potential_files = []
  for fName in os.listdir(path):
    if modelName in fName: #filter for the proper model
      if '.out' in fName:
        if 'tdm.out' not in fName:
          potential_files += [fName]

  dataList = []
  #Now I need to parse all files to get the data
  for fName in potential_files:
    toOpen_indep = open(os.path.join(path, fName))
    try:
      toOpen_dep = open(os.path.join(path, fName.replace(".out", ".tdm.out")))
      dataList += [parser.Parser(toOpen_indep.read(), tdm=toOpen_dep.read())]
    except:
      item = parser.Parser(toOpen_indep.read())
      dataList += [item]
    toOpen_indep.close()

  return dataList