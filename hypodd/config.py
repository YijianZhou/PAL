""" Configure file for hypoDD interface
"""
import os
import numpy as np


class Config(object):
  def __init__(self):

    self.ctlg_code = 'example_pal_hyp-ct'
    self.fsta = 'input/example_pal.sta'
    self.fpha = 'input/example_pal_hyp_full.pha'
    self.dep_corr = 5 # avoid air quake
    self.ot_range = '20190704-20190725'
    self.lat_range = [35.4,36.1]
    self.lon_range = [-117.85,-117.25]
    self.num_grids = [1,1] # x,y (lon, lat)
    self.xy_pad = [0.06,0.05] # degree
    self.num_workers = 10
    self.keep_grids = False

