""" Configure file for hypoDD interface
"""
import os
import numpy as np

class Config(object):
  def __init__(self):

    # 1. format input
    self.fsta_in = 'input/example.sta'
    self.fsta_out = 'input/station.dat'
    self.fpha_in = 'input/example_pad_hyp_all.pha'
    self.fpha_out = 'input/phase.dat'
    self.dep_corr = 5 # avoid air quake

    # 2. format output
    self.out_ctlg = 'output/example_pad_hyp-ct.ctlg'
    self.out_pha = 'output/example_pad_hyp-ct.pha'
    self.out_pha_all = 'output/example_pad_hyp-ct_all.pha'

