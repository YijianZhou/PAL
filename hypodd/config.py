""" Configure file for hypoDD interface
"""
import os
import numpy as np

class Config(object):
  def __init__(self):

    # 1. mk_sta
    self.fsta_in = 'input/sm_station.csv'
    self.fsta_out = 'input/station.dat'

    # 2. mk_phs
    self.fpha_in = 'input/sm_ai_hyp.pha'
    self.fpha_out = 'input/phase.dat'
    self.dep_corr = 0. # avoid air quake

    # 3. reloc2csv
    self.out_ctlg = 'output/sm_ai_reloc.ctlg'
    self.out_pha = 'output/sm_ai_reloc.pha'

