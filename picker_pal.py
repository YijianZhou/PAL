import numpy as np
from scipy.stats import kurtosis

class STA_LTA_Kurtosis(object):
  """ STA/LTA & kurtosis-based P&S Picker
    trigger picker: Z chn STA/LTA reach trig_thres
    --> pick P: find STA/LTA peak within p_win
    --> pick S: find kurtosis peak winthin s_win
  Inputs
    stream: obspy.stream obj (3 chn, [e, n, z])
    win_sta, win_lta: win for sta/lta (det, p, s)
    trig_thres: threshold to trig picker
    p_win, s_win: win len for searching P & S
    pca_win: time win for calc pca filter
    pca_range: time range for pca filter
    win_kurt: win for calc kurtosis
    snr_ratio_thres: max value of snr ratio after peak rm
    amp_ratio_thres: max value of amp ratio of P/P_tail & S
    fd_trhes: min value of dominant frequency
    amp_win: time win to get S amplitude
    det_gap: time gap between detections
    to_prep: whether preprocess stream
    freq_band: frequency band for phase picking
    *note: all time-related params are in sec
  Outputs
    output to file or picks (struct np.array)
  Usage
    import picker_pal
    picker = picker_pal.STA_LTA_Kurtosis()
    picks = picker.pick(stream)
  """
  def __init__(self, 
               win_sta         = [.8, 0.4, 1.],
               win_lta         = [6., 2., 2.],
               trig_thres      = 12.,
               p_win           = [.5, 1.],
               s_win           = 10,
               pca_win         = 1.,
               pca_range       = [0., 2],
               win_kurt        = [5.,1.],
               fd_thres        = 2.5, 
               snr_ratio_thres = 10, 
               amp_ratio_thres = [10,2], 
               amp_win         = [1.,5.],
               det_gap         = 5.,
               to_prep         = True,
               freq_band       = [1., 40]):
    self.win_sta = win_sta
    self.win_lta = win_lta
    self.trig_thres = trig_thres
    self.p_win = p_win
    self.s_win = s_win
    self.pca_win = pca_win
    self.pca_range = pca_range
    self.win_kurt = win_kurt
    self.fd_thres = fd_thres
    self.snr_ratio_thres = snr_ratio_thres
    self.amp_ratio_thres = amp_ratio_thres
    self.amp_win = amp_win
    self.det_gap = det_gap
    self.to_prep = to_prep
    self.freq_band = freq_band

  def pick(self, stream, out_file=None):
    # set output format for picks
    dtype = [('net_sta','O'),
             ('sta_ot','O'),
             ('tp','O'),
             ('ts','O'),
             ('s_amp','O')]
    # preprocess & extract data
    if len(stream)!=3: return np.array([], dtype=dtype)
    if self.to_prep: stream = self.preprocess(stream, self.freq_band)
    if len(stream)!=3: return np.array([], dtype=dtype)
    min_npts = min([len(trace) for trace in stream])
    st_data = np.array([trace.data[0:min_npts] for trace in stream])
    # get header
    head = stream[0].stats
    net_sta = '.'.join([head.network, head.station])
    samp_rate = head.sampling_rate
    start_time, end_time = head.starttime, head.endtime
    # sec to points
    win_sta_npts   = [int(samp_rate * win) for win in self.win_sta]
    win_lta_npts   = [int(samp_rate * win) for win in self.win_lta]
    p_win_npts     = [int(samp_rate * win) for win in self.p_win]
    s_win_npts     =  int(samp_rate * self.s_win)
    pca_win_npts   =  int(samp_rate * self.pca_win)
    pca_range_npts = [int(samp_rate * win) for win in self.pca_range]
    win_kurt_npts  = [int(samp_rate * win) for win in self.win_kurt]
    amp_win_npts   = [int(samp_rate * win) for win in self.amp_win]
    det_gap_npts   =  int(samp_rate * self.det_gap)
    # pick P and S
    picks = []
    # 1. trig picker
    print('1. triggering phase picker')
    cf_trig = self.calc_sta_lta(st_data[2]**2, win_lta_npts[0], win_sta_npts[0])
    trig_index = np.where(cf_trig > self.trig_thres)[0]
    slide_idx = 0
    # 2. phase picking
    print('2. picking phase:')
    for _ in trig_index:
        trig_idx = trig_index[slide_idx]
        if trig_idx < p_win_npts[0] + max(win_lta_npts):
            slide_idx += 1; continue
        # 2.1 pick P with STA/LTA
        p_idx0 = trig_idx - p_win_npts[0] - win_lta_npts[1]
        p_idx1 = trig_idx + p_win_npts[1] + win_sta_npts[1]
        data_p = st_data[2,p_idx0:p_idx1]**2
        cf_p = self.calc_sta_lta(data_p, win_lta_npts[1], win_sta_npts[1])
        tp0_idx = np.argmax(cf_p) + p_idx0
        # check potential glitch
        data_clean = st_data[2,p_idx0:p_idx1].copy()
        peak_idx0 = np.argmax(abs(data_clean))
        peak_idx1 = peak_idx0 + self.find_first_peak(data_clean[peak_idx0:])
        peak_idx0 -= self.find_second_peak(data_clean[0:peak_idx0][::-1])
        peak_idx1 += self.find_second_peak(data_clean[peak_idx1:])+1
        try: data_clean[peak_idx0:peak_idx1] = st_data[2, p_idx0+peak_idx1:p_idx0+2*peak_idx1-peak_idx0]
        except: break
        cf_p_clean = self.calc_sta_lta(data_clean**2, win_lta_npts[1], win_sta_npts[1])
        snr_ratio = np.amax(cf_p) / np.amax(cf_p_clean)
        if snr_ratio>self.snr_ratio_thres: 
            # next detected phase
            rest_det = np.where(trig_index > tp0_idx + det_gap_npts)[0] 
            if len(rest_det)==0: break
            slide_idx = rest_det[0]; continue
        # refine initial pick on waveform
        tp_idx = tp0_idx - self.find_second_peak(data_p[0:tp0_idx-p_idx0][::-1])
        tp = start_time + tp_idx/samp_rate
        # 2.2 pick S 
        # 2.2.1 pca for amp_peak
        if len(st_data[0]) < tp_idx + s_win_npts: break
        s_idx0 = tp_idx - pca_range_npts[0]
        s_idx1 = max(tp_idx + s_win_npts, tp_idx + pca_range_npts[1])
        data_s = np.sum(st_data[0:2, s_idx0:s_idx1]**2, axis=0)**0.5
        pca_filter = self.calc_pca_filter(st_data, tp_idx, pca_range_npts, pca_win_npts)
        data_s[0:len(pca_filter)] *= pca_filter
        dt_peak = max(np.argmax(data_s)+1, pca_win_npts+1)
        # 2.2.2 long_win kurt --> t_max
        s_idx0 = tp_idx + dt_peak//2 - win_kurt_npts[0]
        s_idx1 = tp_idx + dt_peak
        data_s = np.sum(st_data[0:2, s_idx0:s_idx1]**2, axis=0)
        data_s /= np.amax(data_s)
        kurt_long = self.calc_kurtosis(data_s, win_kurt_npts[0])
        # 2.2.3 STA/LTA --> t_min
        s_idx0 = tp_idx + dt_peak//2 - win_lta_npts[2]
        s_idx1 = tp_idx + dt_peak + win_sta_npts[2]
        data_s = np.sum(st_data[0:2, s_idx0:s_idx1]**2, axis=0)
        cf_s = self.calc_sta_lta(data_s, win_lta_npts[2], win_sta_npts[2])[win_lta_npts[2]:]
        # 2.2.4 pick S on short_win kurt
        dt_max = np.argmax(kurt_long) # relative to (tp_idx + dt_peak//2)
        dt_max -= self.find_first_peak(kurt_long[0:dt_max+1][::-1])
        dt_min = np.argmax(cf_s) # relative to (tp_idx + dt_peak//2)
        # if kurt_long not stable, use STA/LTA
        if dt_min>=dt_max: 
            ts0_idx = tp_idx + dt_peak//2 + dt_min
            ts_idx = ts0_idx - self.find_second_peak(data_s[0:dt_min+win_lta_npts[2]][::-1])
        # else, pick peak of kurt_short
        else:
            s_idx0 = tp_idx + dt_peak//2 + dt_min - win_kurt_npts[1]
            s_idx1 = tp_idx + dt_peak//2 + dt_max
            data_s = np.sum(st_data[0:2, s_idx0:s_idx1]**2, axis=0)
            data_s /= np.amax(data_s)
            kurt_short = self.calc_kurtosis(data_s, win_kurt_npts[1])
            kurt_max = np.argmax(kurt_short) if np.argmax(kurt_short)>0 else dt_max-dt_min
            ts0_idx = tp_idx + dt_peak//2 + dt_min + kurt_max
            ts_idx = ts0_idx - self.find_second_peak(data_s[0:s_idx0+win_kurt_npts[1]+kurt_max][::-1])
        ts = start_time + ts_idx/samp_rate if ts_idx>tp_idx else start_time + ts0_idx/samp_rate
        # 3. get related S amplitude
        data_amp = st_data[:, tp_idx-amp_win_npts[0] : ts_idx+amp_win_npts[1]].copy()
        s_amp = self.get_s_amp(data_amp, samp_rate)
        # 4. get p_snr
        p_snr = np.amax(cf_trig[p_idx0:p_idx1])
        # 5. calc dominant frequency & amp ratio
        st = stream.slice(tp, max(tp+(ts-tp)/2, tp+self.pca_win)).copy()
        fd = max([self.calc_freq_dom(tr.data, samp_rate) for tr in st])
        A1 = np.array([np.amax(tr.data)-np.amin(tr.data) for tr in stream.slice(tp, tp+(ts-tp)/2)])
        A2 = np.array([np.amax(tr.data)-np.amin(tr.data) for tr in stream.slice(tp+(ts-tp)/2, ts)])
        A3 = np.array([np.amax(tr.data)-np.amin(tr.data) for tr in stream.slice(ts, ts+(ts-tp)/2)])
        A12 = min([A1[ii]/A2[ii] for ii in range(3)])
        A13 = min([A1[ii]/A3[ii] for ii in range(3)])
        # output picks
        if fd>self.fd_thres and A12<self.amp_ratio_thres[0] and A13<self.amp_ratio_thres[1]:
            print('{}, {}, {}'.format(net_sta, tp, ts))
            sta_ot = self.calc_ot(tp, ts)
            picks.append((net_sta, sta_ot, tp, ts, s_amp))
            if out_file: 
                qual_code = '{:.1f},{:.1f},{:.1f},{:.1f},{:.1f}'.format(p_snr, fd, snr_ratio, A12, A13)
                out_file.write('{},{},{},{},{},{}\n'.format(net_sta, sta_ot, tp, ts, s_amp, qual_code))
        # next detected phase
        rest_det = np.where(trig_index > max(trig_idx,ts_idx,tp_idx) + det_gap_npts)[0]
        if len(rest_det)==0: break
        slide_idx = rest_det[0]
    # convert to structed np.array
    return np.array(picks, dtype=dtype)

  # calc STA/LTA for a trace of data (abs or square)
  def calc_sta_lta(self, data, win_lta_npts, win_sta_npts):
    npts = len(data)
    if npts < win_lta_npts + win_sta_npts:
        print('input data too short!')
        return np.zeros(1)
    sta = np.zeros(npts)
    lta = np.ones(npts)
    data_cum = np.cumsum(data)
    sta[:-win_sta_npts] = data_cum[win_sta_npts:] - data_cum[:-win_sta_npts]
    sta /= win_sta_npts
    lta[win_lta_npts:]  = data_cum[win_lta_npts:] - data_cum[:-win_lta_npts]
    lta /= win_lta_npts
    sta_lta = sta/lta
    sta_lta[0:win_lta_npts] = 0.
    sta_lta[np.isinf(sta_lta)] = 0.
    sta_lta[np.isnan(sta_lta)] = 0.
    return sta_lta

  # calc P wave filter
  def calc_pca_filter(self, data, idx_p, pca_range_npts, pca_win_npts):
    p_mat = data[:, idx_p : idx_p + pca_win_npts]
    p_r, p_v = self.calc_pol(p_mat)
    idx_range = range(idx_p - pca_range_npts[0],
                      idx_p + pca_range_npts[1])
    pca_filter = np.zeros(len(idx_range))
    for i, idx in enumerate(idx_range):
        s_mat = data[:, idx : idx + pca_win_npts]
        s_r, s_v = self.calc_pol(s_mat)
        abs_cos = abs(np.dot(p_v, s_v))
        pca_filter[i] = 1 - s_r * abs_cos
    return pca_filter

  # calc pol_rate & pol_vec
  def calc_pol(self, mat):
    cov = np.cov(mat)
    eig_val, eig_vec = np.linalg.eig(cov)
    lam1  = abs(np.amax(eig_val))
    lam23 = abs(np.sum(eig_val) - lam1)
    pol_rate = 1 - (0.5 * lam23 / lam1)
    pol_vec = eig_vec.T[np.argmax(eig_val)]
    return pol_rate, pol_vec

  # calculate origin time
  def calc_ot(self, tp, ts):
    vp, vs = 5.9, 3.4
    dist = (ts-tp) / (1/vs - 1/vp)
    tt_p = dist / vp
    return tp - tt_p

  # get S amplitide
  def get_s_amp(self, velo, samp_rate):
    # remove mean
    velo -= np.reshape(np.mean(velo, axis=1), [velo.shape[0],1])
    # velocity to displacement
    disp = np.cumsum(velo, axis=1)
    disp /= samp_rate
    return np.amax(abs(np.sum(disp**2, axis=0)))**0.5

  # calc dominant frequency
  def calc_freq_dom(self, data, samp_rate):
    npts = len(data)
    if npts//2==0: return 0
    data -= np.mean(data)
    psd = abs(np.fft.fft(data))**2
    psd = psd[:npts//2]
    return np.argmax(psd) * samp_rate / npts

  # calc kurtosis trace
  def calc_kurtosis(self, data, win_kurt_npts):
    npts = len(data) - win_kurt_npts + 1
    kurt = np.zeros(npts)
    for i in range(npts):
        kurt[i] = kurtosis(data[i:i+win_kurt_npts])
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
    if len(neg_idx)==0 or len(pos_idx)==0: return 0
    first_peak = max(neg_idx[0], pos_idx[0])
    neg_peak = neg_idx[neg_idx>first_peak]
    pos_peak = pos_idx[pos_idx>first_peak]
    if len(neg_peak)==0 or len(pos_peak)==0: return first_peak
    return max(neg_peak[0], pos_peak[0])

  def preprocess(self, stream, freq_band, max_gap=5.):
    # time alignment
    start_time = max([trace.stats.starttime for trace in stream])
    end_time = min([trace.stats.endtime for trace in stream])
    if start_time > end_time: return []
    stream = stream.slice(start_time, end_time, nearest_sample=True)
    # remove nan & inf
    for trace in stream:
        trace.data[np.isnan(trace.data)] = 0
        trace.data[np.isinf(trace.data)] = 0
    # fill data gap
    max_gap_npts = int(max_gap*stream[0].stats.sampling_rate)
    for trace in stream:
        data = trace.data
        npts = len(data)
        data_diff = np.diff(data)
        gap_idx = np.where(data_diff==0)[0]
        gap_list = np.split(gap_idx, np.where(np.diff(gap_idx)!=1)[0] + 1)
        gap_list = [gap for gap in gap_list if len(gap)>=3]
        num_gap = len(gap_list)
        for ii,gap in enumerate(gap_list):
            idx0, idx1 = max(0, gap[0]-1), min(npts-1, gap[-1]+1)
            if ii<num_gap-1: idx2 = min(idx1+(idx1-idx0), idx1+max_gap_npts, gap_list[ii+1][0])
            else: idx2 = min(idx1+(idx1-idx0), idx1+max_gap_npts, npts-1)
            if idx1==idx2: continue
            if idx2==idx1+(idx1-idx0): data[idx0:idx1] = data[idx1:idx2]
            else:
                num_tile = int(np.ceil((idx1-idx0)/(idx2-idx1)))
                data[idx0:idx1] = np.tile(data[idx1:idx2], num_tile)[0:idx1-idx0]
        trace.data = data
    # filter
    stream.detrend('demean').detrend('linear').taper(max_percentage=0.05, max_length=5.)
    freq_min, freq_max = freq_band
    if freq_min and freq_max:
        return stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max)
    elif not freq_max and freq_min:
        return stream.filter('highpass', freq=freq_min)
    elif not freq_min and freq_max:
        return stream.filter('lowpass', freq=freq_max)
    else:
        print('filter type not supported!'); return []

