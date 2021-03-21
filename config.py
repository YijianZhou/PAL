""" Configure file
"""
import data_pipeline as dp

class Config(object):
  def __init__(self):

    # 1. picker params
    self.win_sta     = 0.8         # pick win for STA/LTA
    self.win_lta     = 6.          # pick win for STA/LTA
    self.trig_thres  = 12.         # threshold to trig picker (by energy)
    self.p_win       = [.5,1.2]    # win len for P picking
    self.s_win       = 10.         # win len for S picking
    self.pca_win     = 1.          # win len for PCA filter
    self.pca_range   = [0.,2.]     # time range to apply PCA filter
    self.fd_thres    = 2.5         # min value of dominant frequency
    self.amp_win     = [1.,4.]     # time win to get S amplitude
    self.det_gap     = 5.          # time gap between detections
    self.freq_band   = [2,None]    # frequency band for picking
    self.win_kurt    = 1.          # win for kurtosis calc
    self.amp_thres   = [3, 5]      # kurt backward and forward criteria
    self.peak_gap    = 0.1         # sec, for kurt tail

    # 2. assoc params
    self.ot_dev      = 2.          # time deviation for ot assoc
    self.max_res     = 1.5         # max P res for loc assoc
    self.assoc_num   = 4           # num of stations to assoc
    self.edge_width  = 0.2         # ratio of edge width relative to sta range
    self.xy_grid     = 0.02        # lateral grid width, in degree

    # 3. data interface
    self.get_data_dict = dp.get_data_dict
    self.get_sta_dict = dp.get_sta_dict
    self.get_picks = dp.get_picks
    self.read_data = dp.read_data

