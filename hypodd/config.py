"""
configure file for Hypoinverse input maker
"""
import os
import numpy as np

class Config(object):
  def __init__(self):

    net = 'zsy'
    # 1. mk_sta
    self.fsta_in = '/data2/ZSY_SAC/header/station_ZSY.dat'
    self.fsta_out = 'input/station.dat'

    # 2. mk_phs
    self.fpha_in = '../hypoinverse/output/%s.pha'%net
    self.fpha_out = 'input/phase.dat'
    self.dep_corr = 5. # avoid air quake

    # 3. reloc2csv
    self.out_csv = 'output/%s_reloc.csv'%net
    if not os.path.exists('input'): os.makedirs('input')
    if not os.path.exists('output'): os.makedirs('output')
