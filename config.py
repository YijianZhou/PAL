""" Configure file for PAD.
"""
import data_pipeline as dp
import pickers
import associators

class Config(object):
  def __init__(self):

    # 1. picker params
    self.pick_win    = [10., 1.]      # pick win for STA/LTA
    self.trig_thres  = 15.            # threshold to trig picker (by energy)
    self.pick_thres  = 0.96           # threshold for picking
    self.p_win       = [1., 1.]       # win len for P picking
    self.s_win       = [0, 15.]       # win len for S picking
    self.pca_win     = 1.             # # win len for PCA filter
    self.pca_rng     = [0, 2.5]       # time range to apply PCA filter
    self.fd_thres    = 2.5            # min value of dominant frequency
    self.amp_win     = 5.             # time win to get S amplitude
    self.det_gap     = 5.             # time gap between detections
    self.freq_band   = ['highpass',1] # frequency band for ppk

    # 2. assoc params
    self.ot_dev      = 3.             # time deviation for ot assoc
    self.ttp_dev     = 2.             # time deviation for ttp assoc
    self.assoc_num   = 4              # num of stations to assoc

    # 3. loc params
    self.side_width  = 0.2            # ratio of sides relative to sta range
    self.xy_grid     = 0.02           # lateral grid width, in degree
    self.resp_dict   = {'ZSY': 3.02e8,
                         'YN': 1.67785e9,
                        'XLS': 1/1.6e-9,
                        'G70': 1/1.6e-9,
                        'G4Z': 1/1.6e-9,
                        'G4V': 1/1.6e-9,
                        'REF': 1/1.6e-9
                       }     # instrumental response (cnt/m/s)

    # 4. define func
    self.get_data = dp.get_data
    self.get_picks = dp.get_picks
    self.num_proc = 5
    sta_dict = dp.get_sta('/data2/ZSY_SAC/header/station_ZSY.dat')
    self.picker = pickers.Trad_PS(\
                    trig_thres = self.trig_thres,
                    s_win = self.s_win)
    self.associator = associators.TS_Assoc(\
                        sta_dict, self.resp_dict,
                        assoc_num = self.assoc_num,
                        ot_dev = self.ot_dev,
                        ttp_dev = self.ttp_dev)

