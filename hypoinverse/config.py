"""
configure file for Hypoinverse input maker
"""
import numpy as np

class Config(object):
  def __init__(self):

    net = 'xls'
    # 1. mk_sta
    self.fsta_in = '/data3/XLS_SAC/header/station_XLS.dat'
#    self.fsta_in = '/data2/ZSY_SAC/header/station_ZSY.dat'
    self.fsta_out = 'input/xj.sta'
    self.lat_code = 'N'
    self.lon_code = 'E'

    # 2. mk_phs
    self.fpha_in = '../output/%s/%s_201501-201906_pad.pha'%(net, net)
#    self.fpha_in = '../output/%s/%s_201609-201901_pad.pha'%(net, net)
    self.fpha_out = 'input/%s.phs'%net
    self.mag_corr = 2. # hypoInv do not support neg mag

    # 3. sum2csv
    self.ref_ele = 2.7 # ref ele for CRE mod (max sta ele)
    self.grd_ele = 2. # typical station elevation
    self.fsums = 'output/%s_*.sum'%net
    self.out_csv = 'output/%s.csv'%net
    self.out_sum = 'output/%s.sum'%net
    self.out_bad = 'output/%s_bad.csv'%net
    self.out_good = 'output/%s_good.csv'%net
    self.out_pha = 'output/%s.pha'%net

    # 4. run_hyp
    self.ztr_rng = np.arange(0,40,1)
    self.fhyp_temp = 'temp.hyp'
    self.fhyp = 'input/%s.hyp'%net
    self.pmod = 'input/xj_p_wang.cre'
    self.smod = 'input/xj_s_wang.cre'
    self.net = net

