""" Configure file
"""
import data_pipeline as dp
import pickers
import associators

class Config(object):
  def __init__(self):

    # 1. picker params
    self.pick_win    = [10.,1.]       # pick win for STA/LTA
    self.trig_thres  = 15.            # threshold to trig picker (by energy)
    self.pick_thres  = 0.999          # threshold for picking
    self.p_win       = [1.,1.]        # win len for P picking
    self.s_win       = [0.,12.]       # win len for S picking
    self.pca_win     = 1.             # win len for PCA filter
    self.pca_range   = [0,2.5]        # time range to apply PCA filter
    self.fd_thres    = 2.5            # min value of dominant frequency
    self.amp_win     = [1.,5.]        # time win to get S amplitude
    self.det_gap     = 5.             # time gap between detections
    self.freq_band   = ['highpass',1] # frequency band for ppk

    # 2. assoc params
    self.ot_dev      = 2.             # time deviation for ot assoc
    self.max_res     = 1.5            # max P res for loc assoc
    self.assoc_num   = 4              # num of stations to assoc
    self.edge_width  = 0.2            # ratio of edge width relative to sta range
    self.xy_grid     = 0.02           # lateral grid width, in degree

    # 3. data interface
    self.get_data_dict = dp.get_rc_data
    self.get_sta_dict = dp.get_sta_dict
    self.get_picks = dp.get_picks
    self.read_data = dp.read_rc_data
    self.preprocess = dp.preprocess

