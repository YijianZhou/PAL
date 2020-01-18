""" Configure file
"""
import data_pipelines as dps
import pickers
import associators

class Config(object):
  def __init__(self):

    # 1. picker params
    pick_win    = [10., 1.]      # pick win for STA/LTA
    trig_thres  = 15.            # threshold to trig picker (by energy)
    pick_thres  = 0.96           # threshold for picking
    p_win       = [1., 1.]       # win len for P picking
    s_win       = [0, 15.]       # win len for S picking
    pca_win     = 1.             # win len for PCA filter
    pca_rng     = [0, 2.5]       # time range to apply PCA filter
    fd_thres    = 2.5            # min value of dominant frequency
    amp_win     = 5.             # time win to get S amplitude
    det_gap     = 5.             # time gap between detections
    freq_band   = ['highpass',1] # frequency band for ppk

    # 2. assoc params
    ot_dev      = 3.             # time deviation for ot assoc
    ttp_dev     = 2.             # time deviation for ttp assoc
    assoc_num   = 4              # num of stations to assoc

    # 3. loc params
    side_width  = 0.2            # ratio of sides relative to sta range
    xy_grid     = 0.02           # lateral grid width, in degree
    # sta info
    sta_file    = '/data2/ZSY_SAC/header/station_ZSY.dat'
    resp_dict   = {'ZSY': 3.02e8,
                    'YN': 1.67785e9,
                   'XLS': 1/1.6e-9,
                  }              # instrumental gain (cnt/m/s)

    # 4. define func
    dp = dps.Data(resp_dict)
    self.get_data_dict = dp.get_data_dict
    self.get_picks = dp.get_picks
    self.read_data = dp.read_data
    sta_dict = dp.get_sta_dict(sta_file)
    self.picker = pickers.Trad_PS(\
                    trig_thres = trig_thres,
                    s_win = s_win)
    self.associator = associators.TS_Assoc(\
                        sta_dict,
                        assoc_num = assoc_num,
                        ot_dev = ot_dev,
                        ttp_dev = ttp_dev)

