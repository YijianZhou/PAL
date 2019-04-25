import numpy as np

class Trad_PS(object):

  """ Traditional P&S Picker: STA/LTA for P; PCA + STA/LTA for S
  Algorithm
    trigger picker: Z chn STA/LTA reach thres
    --> pick P: find within p_win
    --> pick S: PCA filter & find winthin s_win
  Inputs
    stream: obspy.stream obj (3 chn, [e, n, z])
    pick_win: win len for STA/LTA ([lwin, swin])
    trig_thres: threshold to trig picker
    pick_thres: threshold for picking (0 to 1.)
    p_win: win len for pick detected P
    s_win: win len for S arrivla searching
    pca_win: win len for calc PCA
    s_stride: time stride for S picking
    amp_win: time win to get S amplitude
    *all time related params are in sec
  Outputs
    all picks in the stream, and header info
  Usage
    import pickers
    picker = pickers.Trad_PS()
    picks = picker.pick()
  """

  def __init__(self, 
               pick_win    = [10., 1.],
               trig_thres  = 20., 
               pick_thres  = 0.96,
               p_win       = 2.,
               pca_win     = 1.,
               s_win       = 20.,
               s_stride    = 0.05,
               amp_win     = 5.):

    # change sec to points for time params
    self.pick_win    = [int(100*pick_win[0]), int(100*pick_win[1])] 
    self.idx_shift   = self.pick_win
    self.trig_thres  = trig_thres
    self.pick_thres  = pick_thres
    self.p_win       = [int(100*p_win), int(100*p_win)]
    self.s_win       = [-1*100, int(100*s_win)] #TODO
    self.pca_win     = int(100*pca_win)
    self.s_stride    = int(100*s_stride)
    self.amp_win     = int(100*amp_win)


  def pick(self, stream):

    # set output format
    dtype = [('network','O'),
             ('station','O'),
             ('sta_lon','O'),
             ('sta_lat','O'),
             ('org_t0','O'),
             ('p_arr','O'),
             ('s_arr','O'),
             ('s_amp','O'),
             ('p_snr','O'),
             ('s_snr','O')]

    # time alignment
    start_time = max([trace.stats.starttime for trace in stream])
    end_time   = min([trace.stats.endtime for trace in stream])
    if start_time>end_time: return np.array([],dtype=dtype)
    stream = stream.slice(start_time, end_time)
    if len(stream)!=3: return np.array([],dtype=dtype)

    # read header
    head = stream[0].stats
    net     = head.network
    sta     = head.station
    sta_lon = head.sac.stlo
    sta_lat = head.sac.stla

    # read data
    stream.detrend('demean').detrend('linear').\
           taper(max_percentage=0.05).\
           filter('highpass', freq=1.)
    datax = stream[0].data
    datay = stream[1].data
    dataz = stream[2].data
    data  = [datax, datay, dataz]    

    # pick P and S
    picks = []
    print('-'*40)
    # 1. trig picker
    print('triggering phase picker: {}, {}'.format(sta, start_time))
    cf_trig = self.calc_cf(dataz, self.pick_win)
    self.cf_trig=cf_trig #TODO
    trig_ppk = np.where(cf_trig > self.trig_thres)[0]
    slide_idx = 0
    print('picking phase:')
    for _ in trig_ppk:

        # 2. pick detected P
        idx_trig = trig_ppk[slide_idx]
        if idx_trig -self.p_win[0] -self.idx_shift[0] <0:
            slide_idx+=1; continue
        data_p = dataz[idx_trig -self.p_win[0] -self.idx_shift[0]
                      :idx_trig +self.p_win[1] +self.idx_shift[1]]
        cf_p = self.calc_cf(data_p, self.pick_win)
        idx_p = idx_trig -self.idx_shift[0] -self.p_win[0] +\
                np.where(cf_p >= self.pick_thres *np.amax(cf_p))[0][0]
        tp = start_time + idx_p/100

        # 3. pick related S
        if min(len(data[0]), len(data[1]), len(data[2]))\
           < idx_p +self.s_win[1] +self.pca_win: break
        data_s = self.calc_filter(data, idx_p, self.s_stride)
        cf_s  = self.calc_cf(data_s,  self.pick_win, decim=self.s_stride)
        self.cf_s = cf_s #TODO
        self.data_s = data_s #TODO

        # trig S picker and pick
        s_trig = np.argmax(data_s)
        s_rng0 = min(s_trig, int(s_trig +self.idx_shift[0]/self.s_stride)//2)
        s_rng1 = max(s_trig, int(s_trig +self.idx_shift[0]/self.s_stride)//2)
        if s_rng0==s_rng1: s_rng1+=1
        cf_s = cf_s[s_rng0:s_rng1]
        idx_s = idx_p -self.idx_shift[0] -self.s_win[0]\
                +self.s_stride*(s_rng0+np.argmax(cf_s))
        ts = start_time + idx_s/100

        # get related S amplitude
        ampx = self.get_amp(datax[idx_s :idx_s+self.amp_win])
        ampy = self.get_amp(datay[idx_s :idx_s+self.amp_win])
        ampz = self.get_amp(dataz[idx_s :idx_s+self.amp_win])
        amp = np.sqrt(ampx**2 + ampy**2 + ampz**2)

        # get p_anr and s_anr
        p_snr = np.amax(cf_p)
        s_snr = np.amax(cf_s)

        # output
        print('{}, {}, {}'.format(sta, tp, ts))
        if not tp>ts:
            ot0 = self.est_ot(tp, ts) # est. of ot for assoc
            picks.append((net, sta, sta_lon, sta_lat, ot0, tp, ts, amp, p_snr, s_snr))

        # next detected phase
        rest_det = np.where(trig_ppk >\
                   max(idx_trig, idx_s+5*100, idx_p+5*100))[0] #TODO
        if len(rest_det)==0: break
        slide_idx = rest_det[0]

    # convert to structed np.array
    return np.array(picks, dtype=dtype)


  def calc_cf(self, data, win, decim=1.):
    """  calc character func of STA/LTA for a single trace
    Inputs
        data: input trace data, in np.array
        win: win len for STA/LTA, [lwin, swin], in points
        decim: decimate rate
    Outputs
        cf: character function
    """
    lwin = int(win[0] /decim)
    swin = int(win[1] /decim)
    if len(data)<lwin+swin:
        print('input data too short!')
        return np.zeros(1)
    sta = np.zeros(len(data))
    lta = np.ones(len(data))
    # use energy
    data = np.cumsum(data**2)
    # Compute the STA and the LTA
    sta[:-swin] = data[swin:] - data[:-swin]
    sta /= swin
    lta[lwin:]  = data[lwin:] - data[:-lwin]
    lta /= lwin
    # Pad zeros (same out size as data)
    sta[:lwin] = 0
    cf = sta/lta
    # avoid bad points
    cf[np.isinf(cf)] = 0.
    cf[np.isnan(cf)] = 0.
    return cf


  def calc_filter(self, data, idx_p, stride):
    """ calc S filter by PCA
    Inputs:
        data: input 3 chn array
        idx_p: idx for P in data
        stride: stride for the filter points
    Outputs:
        data_s0, data_s1: filtered S traces (E&N)
    """
    # calc P pol
    p_mat = np.array([data[0][idx_p: idx_p+self.pca_win],
                      data[1][idx_p: idx_p+self.pca_win],
                      data[2][idx_p: idx_p+self.pca_win]])
    p_r, p_evec = self.calc_pol(p_mat)
    # calc S filter
    idx_rng = range(idx_p -self.s_win[0] -self.idx_shift[0], 
                    idx_p +self.s_win[1] +self.idx_shift[1], stride)
    flt = np.zeros(len(idx_rng))
    for i, idx in enumerate(idx_rng):
        s_mat = np.array([data[0][idx: idx+self.pca_win],
                          data[1][idx: idx+self.pca_win],
                          data[2][idx: idx+self.pca_win]])
        s_r, s_evec = self.calc_pol(s_mat)
        u11 = abs(np.dot(p_evec, s_evec))
        flt[i] = s_r *(1-u11)
    data_s = flt* np.sqrt(data[0][idx_rng]**2 + data[1][idx_rng]**2)
    return data_s

  def calc_pol(self, mat):
    """ calc polarization by PCA
    Inputs
        mat: 3chn time win (matrix)
    Outputs
        r: polirization degree
        vec: dominant eig-vector
    """
    cov = np.cov(mat)
    eval, evec = np.linalg.eig(cov)
    # calc pol degree
    lam1  = np.amax(eval)
    lam23 = np.sum(eval) - lam1
    r = 1 - (0.5*lam23/lam1)
    # calc dom vec
    vec = evec.T[np.argmax(eval)]
    return r, vec


  # estimate original time
  def est_ot(self, tp, ts):
    vp = 5.9
    vs = 3.4
    r = (ts-tp) /(1/vs - 1/vp)
    Tp = r/vp
    return tp-Tp

  # get S amplitide
  def get_amp(self, velo):
    # velocity to displacement
    disp = np.zeros(len(velo))
    for i in range(len(velo)-1):
        disp[i+1] = np.sum(velo[0:i])
    disp = 0.01*disp # 100 sps
    return (np.amax(disp) - np.amin(disp)) /2

