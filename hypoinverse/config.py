""" Configure file for HypoInverse interface
"""
import numpy as np

class Config(object):
  def __init__(self):

    self.ctlg_code = 'eg_pal_hyp'
    # 1. mk_sta: format station file
    self.fsta = 'input/example_pal.sta'
    self.lat_code = 'N'
    self.lon_code = 'W'
    # 2. mk_pha: format phase file
    self.fpha = 'input/eg_pal.pha'
    self.mag_corr = 2. # hypoInv do not support neg mag
    # 3. sum2csv: format output files
    self.ref_ele = 2.5 # ref ele for CRE mod (max sta ele)
    self.grd_ele = 1.5 # typical station elevation
    self.fsums = 'output/%s-*.sum'%self.ctlg_code
    self.out_ctlg = 'output/%s.ctlg'%self.ctlg_code
    self.out_pha = 'output/%s.pha'%self.ctlg_code
    self.out_pha_full = 'output/%s_full.pha'%self.ctlg_code
    self.out_sum = 'output/%s.sum'%self.ctlg_code
    self.out_bad = 'output/%s_bad.csv'%self.ctlg_code
    self.out_good = 'output/%s_good.csv'%self.ctlg_code
    # 4. run_hyp
    self.num_workers = 10
    self.ztr_rng = np.arange(0,20,1)
    self.p_wht = 0 # weight code
    self.s_wht = 1
    self.rms_wht = '4 0.3 1 3'
    self.dist_init = '1 50 1 2'
    self.dist_wht = '4 20 1 3'
    self.wht_code = '1 0.6 0.3 0.2'
    self.fhyp_temp = 'temp_hyp/temp_vp-pos.hyp'
    self.pmod = 'input/example_p.cre'
    self.smod = [None, 'input/example_s.cre'][0]
    self.pos = 1.73 # provide smod or pos
    self.get_prt = False
    self.get_arc = False
    self.keep_fsums = False

