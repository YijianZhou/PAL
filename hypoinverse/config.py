""" Configure file for HypoInverse interface
"""
import numpy as np

class Config(object):
  def __init__(self):

    # 1. mk_sta: format station file
    self.fsta_in = 'input/sm_station.csv'
    self.fsta_out = 'input/sm.sta'
    self.lat_code = 'N'
    self.lon_code = 'E'

    # 2. mk_phs: format phase file
    self.fpha_in = 'input/phase_ai.csv'
    self.fpha_out = 'input/sm_ai.phs'
    self.mag_corr = 2. # hypoInv do not support neg mag

    # 3. sum2csv: format output files
    self.ref_ele = 3. # ref ele for CRE mod (max sta ele)
    self.grd_ele = 2. # typical station elevation
    self.fsums = 'output/sm_ai-*.sum'
    self.out_ctlg = 'output/sm_ai_hyp.ctlg'
    self.out_pha = 'output/sm_ai_hyp.pha'
    self.out_sum = 'output/sm_ai.sum'
    self.out_bad = 'output/sm_ai_bad.csv'
    self.out_good = 'output/sm_ai_good.csv'

    # 4. run_hyp
    self.num_workers = 10
    self.ztr_rng = np.arange(0,30,1)
    self.p_wht = 0 # weight code
    self.s_wht = 1
    self.rms_wht = '4 0.3 1 2'
    self.dist_init = '1 40 2 4'
    self.dist_wht = '4 50 1 2'
    self.wht_code = '1 0.5 0.5 0.2'
    self.fhyp_temp = 'temp_hyp/temp_vp-pos.hyp'
    self.pmod = 'input/sm_p_wang.cre'
    self.pos = 1.83
    self.smod = None
    self.get_prt = False
    self.get_arc = False
    self.ctlg_code = 'sm_ai'
    self.keep_fsums = False

