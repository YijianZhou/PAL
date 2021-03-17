import numpy as np
from scipy.stats import kurtosis


class STA_LTA_PCA(object):

  """ STA/LTA based P&S Picker
    trigger picker: Z chn STA/LTA reach trig_thres
    --> pick P: find kurtosis peak within p_win
    --> pick S: PCA filter & find kurtosis peak winthin s_win
  Inputs
    stream: obspy.stream obj (3 chn, [e, n, z])
    pick_win: win len for STA/LTA ([lwin, swin])
    trig_thres: threshold to trig picker
    p_win: win len for pick detected P
    s_win: win len for S arrivla searching
    pca_win: time win for calc pca filter
    pca_range: time range for pca filter
    fd_trhes: minimum value of dominant frequency
    amp_win: time win to get S amplitude
    det_gap: time gap between detections
    to_prep: whether preprocess stream
    freq_band: frequency band for phase picking
    *note: all time-related params are in sec
  Outputs
    all picks in the stream, and header info
  Usage
    import picker_pal
    picker = picker_pal.STA_LTA_PCA()
    picks = picker.pick(stream)
  """

  def __init__(self, 
               win_sta    = 1.,
               win_lta    = 10.,
               trig_thres = 15.,
               p_win      = [.5, 1.2],
               s_win      = 12,
               pca_win    = 1.,
               pca_range  = [0., 2],
               fd_thres   = 2.5,
               amp_win    = [1.,5.],
               win_kurt   = 1.,
               amp_thres  = [3, 5],  
               peak_gap   = 0.1, 
               det_gap    = 5.,
               to_prep    = True,
               freq_band  = [1., None]):

    self.win_sta    = win_sta
    self.win_lta    = win_lta
    self.trig_thres = trig_thres
    self.p_win      = p_win
    self.s_win      = s_win
    self.pca_win    = pca_win
    self.pca_range  = pca_range
    self.win_kurt   = win_kurt
    self.peak_gap   = peak_gap
    self.amp_thres  = amp_thres
    self.fd_thres   = fd_thres
    self.amp_win    = amp_win
    self.det_gap    = det_gap
    self.to_prep    = to_prep
    self.freq_band  = freq_band


  def pick(self, stream, out_file=None):
    # set output format
    dtype = [('net_sta','O'),
             ('sta_ot','O'),
             ('tp','O'),
             ('ts','O'),
             ('s_amp','O'),
             ('p_snr','O'),
             ('freq_dom','O')]
    # preprocess & extract data
    if len(stream)!=3: return np.array([], dtype=dtype)
    if self.to_prep: stream = self.preprocess(stream, self.freq_band)
    if len(stream)!=3: return np.array([], dtype=dtype)
    min_npts = min([len(trace) for trace in stream])
    st_data = np.array([trace.data[0:min_npts] for trace in stream])
    # get header
    head = stream[0].stats
    net_sta = '.'.join([head.network, head.station])
    self.samp_rate = head.sampling_rate
    start_time, end_time = head.starttime, head.endtime
    # sec to points
    self.win_sta_npts   =  int(self.samp_rate * self.win_sta)
    self.win_lta_npts   =  int(self.samp_rate * self.win_lta)
    self.p_win_npts     = [int(self.samp_rate * win) for win in self.p_win]
    self.s_win_npts     =  int(self.samp_rate * self.s_win)
    self.pca_win_npts   =  int(self.samp_rate * self.pca_win)
    self.pca_range_npts = [int(self.samp_rate * win) for win in self.pca_range]
    self.win_kurt_npts  =  int(self.samp_rate * self.win_kurt)
    amp_win_npts        = [int(self.samp_rate * win) for win in self.amp_win]
    det_gap_npts        =  int(self.samp_rate * self.det_gap)
    self.peak_gap_npts  =  int(self.samp_rate * self.peak_gap)

    # pick P and S
    picks = []
    # 1. trig picker
    print('1. triggering phase picker')
    cf_trig = self.calc_sta_lta(st_data[2])
    trig_index = np.where(cf_trig > self.trig_thres)[0]
    slide_idx = 0
    # 2. phase picking
    print('2. picking phase:')
    for _ in trig_index:
        trig_idx = trig_index[slide_idx]
        if trig_idx < self.p_win_npts[0] + self.win_lta_npts: 
            slide_idx += 1; continue
        # 2.1 pick P
        p_idx0 = trig_idx - self.p_win_npts[0] - self.win_kurt_npts
        p_idx1 = trig_idx + self.p_win_npts[1]
        data_p = np.sum(st_data[:,p_idx0:p_idx1]**2, axis=0)
        data_p /=  np.amax(data_p)
        p_idx = self.pick_kurtosis(data_p) + p_idx0
        tp = start_time + p_idx / self.samp_rate
        # 2.2 pick S between tp and tp+s_win
        # pca for s_peak
        if len(st_data[0]) < p_idx + self.s_win_npts: break
        s_idx0 = p_idx - self.pca_range_npts[0]
        s_idx1 = max(p_idx + self.s_win_npts, p_idx + self.pca_range_npts[1])
        data_s = np.sum(st_data[0:2, s_idx0:s_idx1]**2, axis=0)
        pca_filter = self.calc_pca_filter(st_data, p_idx)
        data_s[0:len(pca_filter)] *= pca_filter
        dt_peak = max(np.argmax(data_s)+1, self.pca_win_npts+1)
        # pick s between tp+(s_peak-tp)/2 and s_peak
        s_idx0 = p_idx + dt_peak//2 - self.win_kurt_npts
        s_idx1 = p_idx + dt_peak
        data_s = np.sum(st_data[0:2, s_idx0:s_idx1]**2, axis=0)
        data_s /= np.amax(data_s)
        s_idx = self.pick_kurtosis(data_s) + s_idx0
        ts = start_time + s_idx / self.samp_rate
        # 2.3 get related S amplitude
        s_amp = self.get_s_amp(st_data[:, p_idx-amp_win_npts[0] : s_idx+amp_win_npts[1]])
        # 2.4 get p_anr
        p_snr = np.amax(cf_trig[p_idx0:p_idx1])
        # 2.5 calc dominant frequency
        t0 = min(tp, ts)
        t1 = max(tp+(ts-tp)/2, tp+self.win_sta)
        st = stream.slice(t0,t1)
        fd = max([self.calc_freq_dom(tr.data) for tr in st])
        # output
        print('{}, {}, {}'.format(net_sta, tp, ts))
        if tp<ts and fd>self.fd_thres:
            sta_ot = self.calc_ot(tp, ts)
            picks.append((net_sta, sta_ot, tp, ts, s_amp, p_snr, fd))
            if out_file: 
                pick_line = '{},{},{},{},{},{:.2f},{:.2f}\n'\
                    .format(net_sta, sta_ot, tp, ts, s_amp, p_snr, fd)
                out_file.write(pick_line)
        # next detected phase
        rest_det = np.where(trig_index > max(trig_idx,s_idx,p_idx) + det_gap_npts)[0]
        if len(rest_det)==0: break
        slide_idx = rest_det[0]
    # convert to structed np.array
    return np.array(picks, dtype=dtype)


  def pick_kurtosis(self, data):
    # 1. calc kurtosis 
    kurt = self.calc_kurtosis(data)
    if len(kurt)<self.peak_gap_npts: return len(data)
    # 2. choose the correct kurt peak
    peak_idx = np.argmax(data[self.win_kurt_npts:]) # find kurt_peak before amp_peak
    if peak_idx<self.peak_gap_npts: return self.win_kurt_npts + self.peak_gap_npts
    k_idx = np.argmax(kurt[0:peak_idx+1]) # idx in kurt domain
    amp = self.get_peak_amp(data[self.win_kurt_npts:], k_idx)
    dt1_idx = self.find_first_peak(kurt[0:k_idx][::-1])
    dt2_idx = max(self.find_first_peak(kurt[k_idx:]), self.peak_gap_npts)
    k1_min, k1_max = 0, max(k_idx-dt1_idx, 1)
    k2_min = min(k_idx+dt2_idx, peak_idx-1)
    k2_max = max(k2_min+1, peak_idx+1)
    k1_idx = np.argmax(kurt[k1_min:k1_max]) 
    k2_idx = np.argmax(kurt[k2_min:k2_max]) + k2_min 
    amp1 = self.get_peak_amp(data[self.win_kurt_npts:], k1_idx)
    amp2 = self.get_peak_amp(data[self.win_kurt_npts:], k2_idx)
    if amp2/amp>self.amp_thres[1]: k_idx = k2_idx
    elif amp/amp1<=self.amp_thres[0]: k_idx = k1_idx
    t_data = self.win_kurt_npts + k_idx # idx in data domain
    dt_data = self.find_second_peak(data[0:t_data][::-1])
    return t_data - dt_data


  # calc STA/LTA for a trace of data
  def calc_sta_lta(self, data):
    npts = len(data)
    if npts < self.win_lta_npts + self.win_sta_npts:
        print('input data too short!')
        return np.zeros(1)
    sta = np.zeros(npts)
    lta = np.ones(npts)
    # use energy
    data = np.cumsum(data**2)
    # calc STA and LTA
    sta[:-self.win_sta_npts] = data[self.win_sta_npts:] - data[:-self.win_sta_npts]
    sta /= self.win_sta_npts
    lta[self.win_lta_npts:]  = data[self.win_lta_npts:] - data[:-self.win_lta_npts]
    lta /= self.win_lta_npts
    # pad zeros (same out size as data)
    sta[:self.win_lta_npts] = 0
    sta_lta = sta/lta
    # avoid bad points
    sta_lta[np.isinf(sta_lta)] = 0.
    sta_lta[np.isnan(sta_lta)] = 0.
    return sta_lta


  # calc P wave filter
  def calc_pca_filter(self, data, idx_p):
    p_mat = data[:, idx_p : idx_p + self.pca_win_npts]
    p_r, p_v = self.calc_pol(p_mat)
    # calc filter
    idx_range = range(idx_p - self.pca_range_npts[0],
                      idx_p + self.pca_range_npts[1])
    pca_filter = np.zeros(len(idx_range))
    for i, idx in enumerate(idx_range):
        s_mat = data[:, idx : idx + self.pca_win_npts]
        s_r, s_v = self.calc_pol(s_mat)
        abs_cos = abs(np.dot(p_v, s_v))
        pca_filter[i] = 1 - s_r * abs_cos
    return pca_filter


  # calc pol_rate & pol_vec
  def calc_pol(self, mat):
    cov = np.cov(mat / np.amax(abs(mat)))
    eig_val, eig_vec = np.linalg.eig(cov)
    # calc pol degree
    lam1  = abs(np.amax(eig_val))
    lam23 = abs(np.sum(eig_val) - lam1)
    pol_rate = 1 - (0.5 * lam23 / lam1)
    # get pol vec
    pol_vec = eig_vec.T[np.argmax(eig_val)]
    return pol_rate, pol_vec


  # calculate origin time
  def calc_ot(self, tp, ts):
    vp, vs = 5.9, 3.4
    dist = (ts-tp) / (1/vs - 1/vp)
    tt_p = dist / vp
    return tp - tt_p


  # get S amplitide
  def get_s_amp(self, velo):
    # remove mean
    velo -= np.reshape(np.mean(velo, axis=1), [velo.shape[0],1])
    # velocity to displacement
    disp = np.cumsum(velo, axis=1)
    disp /= self.samp_rate
    return np.amax(np.sum(disp**2, axis=0))**0.5


  # calc dominant frequency
  def calc_freq_dom(self, data):
    npts = len(data)
    if npts//2==0: return 0
    data -= np.mean(data)
    psd = abs(np.fft.fft(data))**2
    psd = psd[:npts//2]
    return np.argmax(psd) * self.samp_rate / npts


  # calc kurtosis trace
  def calc_kurtosis(self, data):
    npts = len(data) - self.win_kurt_npts + 1
    kurt = np.zeros(npts)
    for i in range(npts):
        kurt[i] = kurtosis(data[i:i+self.win_kurt_npts])
    return kurt


  def find_first_peak(self, data):
    npts = len(data)
    if npts<2: return 0
    delta_d = data[1:npts] - data[0:npts-1]
    if min(delta_d)>=0 or max(delta_d)<=0: return 0
    neg_idx = np.where(delta_d<0)[0]
    pos_idx = np.where(delta_d>=0)[0]
    return max(neg_idx[0], pos_idx[0])


  def find_second_peak(self, data):
    npts = len(data)
    if npts<2: return 0
    delta_d = data[1:npts] - data[0:npts-1]
    if min(delta_d)>=0 or max(delta_d)<=0: return 0
    neg_idx = np.where(delta_d<0)[0]
    pos_idx = np.where(delta_d>=0)[0]
    first_peak = max(neg_idx[0], pos_idx[0])
    neg_peak = neg_idx[neg_idx>first_peak]
    pos_peak = pos_idx[pos_idx>first_peak]
    if len(neg_peak)==0 or len(pos_peak)==0: return first_peak
    return max(neg_peak[0], pos_peak[0])


  def get_peak_amp(self, data, t0):
    idx0 = max(t0 - self.find_second_peak(data[0:t0][::-1]), 0)
    idx1 = min(t0 + self.find_second_peak(data[t0:]), len(data)-1)
    if idx0==idx1: return data[idx0]
    else: return np.amax(data[idx0:idx1])


  def preprocess(self, stream, freq_band):
    # time alignment
    start_time = max([trace.stats.starttime for trace in stream])
    end_time = min([trace.stats.endtime for trace in stream])
    if start_time > end_time: return []
    stream = stream.slice(start_time, end_time, nearest_sample=True)
    # filter
    stream.detrend('demean').detrend('linear').taper(max_percentage=0.05, max_length=10.)
    freq_min, freq_max = freq_band
    if freq_min and freq_max:
        return stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
    elif not freq_max and freq_min:
        return stream.filter('highpass', freq=freq_min)
    elif not freq_min and freq_max:
        return stream.filter('lowpass', freq=freq_max)
    else:
        print('filter type not supported!'); return []

