import numpy as np

class TS_Assoc(object):

  """ Associate picks by searching ot (time, T) and ttp (space, S) cluster
  Inputs
    sta_dict: station location dict
    side_width: ratio of sides relative to the station range
    xy_grid: grid width for x-y axis (in degree)
    vp, vs: P&S velo of the uniform model
    ot_dev: time dev for ot assoc
    ttp_dev: time dev for ttp assoc
    assoc_num: min num to alert a detection
    *note: lateral spatial range (x-y) in degree; depth in km; elevation in m
  Usage
    import associators
    associator = associators.TS_Assoc(sta_dict)
    associator.associate(picks, out_ctlg, out_pha)
  """
  
  def __init__(self,
               sta_dict,
               side_width = 0.2,
               xy_grid    = 0.02,
               vp         = 5.9,
               vs         = 3.4,
               ot_dev     = 3.,
               ttp_dev    = 2.,
               assoc_num  = 4):

    self.sta_dict   = sta_dict
    self.side_width = side_width
    self.xy_grid    = xy_grid
    self.vp         = vp
    self.vs         = vs
    self.ot_dev     = ot_dev
    self.ttp_dev    = ttp_dev
    self.assoc_num  = assoc_num
    self.time_table = self.calc_tt()


  def associate(self, picks, out_ctlg=None, out_pha=None):

    # 1. temporal assoc: picks --> event_picks
    event_picks = self.pick2event(picks)

    # 2. spatial assoc: event_pick --> event_loc
    for event_pick in event_picks:
        event_loc, event_pick = self.locate(event_pick)
        if len(event_loc)==0: continue
        # 3. estimate magnitude
        event_loc_mag = self.calc_mag(event_pick, event_loc)
        # write catalog and phase
        self.write_catalog(event_loc_mag, out_ctlg)
        self.write_phase(event_loc_mag, event_pick, out_pha)


  # 1. temporal assoc by ot cluster: picks --> event_picks
  def pick2event(self, picks):
    event_picks = []
    num_picks = len(picks)
    if num_picks==0: print('no event detected'); return event_picks
    picks = np.sort(picks, order='sta_ot')

    # calc num_nbr
    ots = picks['sta_ot']
    num_nbr = np.zeros(num_picks)
    for i, oti in enumerate(ots):
        is_nbr = abs(ots-oti) < self.ot_dev
        num_nbr[i] = sum(is_nbr.astype(float))

    # assoc each cluster
    for _ in range(num_picks):
        if np.amax(num_nbr) < self.assoc_num: break
        # ot assoc
        oti = ots[np.argmax(num_nbr)]
        to_assoc = abs(ots-oti) < self.ot_dev
        event_picks.append(picks[to_assoc])
        num_nbr[to_assoc] = 0
        # renew num_nbr
        to_renew  = abs(ots-oti) >   self.ot_dev
        to_renew *= abs(ots-oti) < 2*self.ot_dev
        ots_todel = ots[to_assoc]
        nbr_todel = [sum(abs(ots_todel-oti) < self.ot_dev) for oti in ots[to_renew]]
        num_nbr[to_renew] -= nbr_todel
    print('associated {} events'.format(len(event_picks)))
    return event_picks


  # 2. spatial assoc by locate event (sta res): event_pick --> event_loc
  def locate(self, event_pick):
    res_ttp = 0 # P travel time res
    tp_num = 0 # true positive
    tp_dict = {}
    num_sta = len(event_pick)
    ot = event_pick['sta_ot'][num_sta//2]
    for pick in event_pick:
        net = pick['net']
        sta = pick['sta']
        net_sta = '.'.join([net, sta])
        ttp_obs = pick['p_arr'] - ot # to travel time
        ttp_pred = self.time_table[net_sta]
        resi = abs(ttp_pred - ttp_obs)
        resi[resi > self.ttp_dev] = 0.
        res_ttp += resi
        is_tp = (resi>0.).astype(float)
        tp_dict[net_sta] = is_tp
        tp_num += is_tp

    # find min res in the map
    num_sta = np.amax(tp_num)
    if num_sta < self.assoc_num: return [],[]
    res_ttp /= tp_num
    res_ttp[tp_num != num_sta] = np.inf
    min_res = np.amin(res_ttp)
    x, y = np.unravel_index(np.argmin(res_ttp), res_ttp.shape)
    lon = self.lon_rng[0] + x * self.xy_grid
    lat = self.lat_rng[0] + y * self.xy_grid

    # find sta phase
    event_pick = [pick for pick in event_pick \
               if tp_dict[pick['net'] +'.'+ pick['sta']][x][y] == 1.]

    # output
    print('locate event: ot {}, lon {:.2f}, lat {:.2f}, res {:.1f}'\
      .format(ot, lon, lat, min_res))
    event_loc = {'evt_ot' : ot, 
                 'evt_lon': round(lon,2), 
                 'evt_lat': round(lat,2),
                 'res': round(min_res,1)}
    return event_loc, event_pick


  # calc time table (init assoc)
  def calc_tt(self):
    # set up
    print('making time table')
    time_table = {}
    dist_grid = 111 * self.xy_grid # in km

    # get x-y range: sta range + side width
    lon, lat = [], []
    for sta_loc in self.sta_dict.values():
        lon.append(sta_loc['sta_lon'])
        lat.append(sta_loc['sta_lat'])
    lon, lat = np.array(lon), np.array(lat)
    # calc side width
    lon_side = self.side_width * (np.amax(lon) - np.amin(lon))
    lat_side = self.side_width * (np.amax(lat) - np.amin(lat))
    # get range
    lon_rng = [np.amin(lon) - lon_side, np.amax(lon) + lon_side]
    lat_rng = [np.amin(lat) - lat_side, np.amax(lat) + lat_side]
    # set x-y axis
    x_rng = range(int((lon_rng[1] - lon_rng[0]) / self.xy_grid))
    y_rng = range(int((lat_rng[1] - lat_rng[0]) / self.xy_grid))

    # calc time table
    for net_sta, sta_loc in self.sta_dict.items():
        # convert to x-y axis
        lon = sta_loc['sta_lon']
        lat = sta_loc['sta_lat']
        ele = sta_loc['sta_ele'] / 1000. # m to km
        sta_x = int((lon - lon_rng[0]) / self.xy_grid)
        sta_y = int((lat - lat_rng[0]) / self.xy_grid)
        # calc P travel time
        ttp = -np.ones([len(x_rng), len(y_rng)])
        for i,x in enumerate(x_rng):
          for j,y in enumerate(y_rng):
            dist = dist_grid * np.sqrt((x-sta_x)**2 + (y-sta_y)**2)
            dep  = ele + 5 # init event depth 5km
            ttp[i,j] = np.sqrt(dep**2 + dist**2) / self.vp
        time_table[net_sta] = ttp
    self.lon_rng = lon_rng
    self.lat_rng = lat_rng
    return time_table


  # calc mag with picks (s_amp)
  def calc_mag(self, event_pick, event_loc):
    num_sta = len(event_pick)
    mag = np.zeros(num_sta)
    dist = np.zeros(num_sta)
    for i,pick in enumerate(event_pick):
        # get sta_loc
        net_sta = '.'.join([pick['net'], pick['sta']])
        sta_loc = self.sta_dict[net_sta]
        # get S amp
        amp = pick['s_amp'] * 1e6 # m to miu m
        # calc epi dist
        dist_lon = 111*(sta_loc['sta_lon'] - event_loc['evt_lon'])
        dist_lat = 111*(sta_loc['sta_lat'] - event_loc['evt_lat'])
        dist[i] = np.sqrt(dist_lon**2 + dist_lat**2) # in km
        mag[i] = np.log10(amp) + np.log10(dist[i])
    event_loc['mag'] = round(np.median(mag),2)
    print('estimated magnitude: {:.1f} delta {:.1f}'\
      .format(np.median(mag), np.std(mag)))
    return event_loc


  # write event loc into catalog
  def write_catalog(self, event_loc, out_ctlg):
    ot  = event_loc['evt_ot']
    lon = event_loc['evt_lon']
    lat = event_loc['evt_lat']
    mag = event_loc['mag']
    res = event_loc['res']
    out_ctlg.write('{},{},{},{},{}\n'.format(ot, lat, lon, mag, res))


  # write sta phase into phase file
  def write_phase(self, event_loc, event_pick, out_pha):
    ot  = event_loc['evt_ot']
    lon = event_loc['evt_lon']
    lat = event_loc['evt_lat']
    mag = event_loc['mag']
    res = event_loc['res']
    out_pha.write('{},{},{},{},{}\n'.format(ot, lat, lon, mag, res))
    for pick in event_pick:
        net   = pick['net']
        sta   = pick['sta']
        p_arr = pick['p_arr']
        s_arr = pick['s_arr']
        s_amp = pick['s_amp']
        p_snr = pick['p_snr']
        s_snr = pick['s_snr']
        out_pha.write('{},{},{},{},{},{:.1f},{:.1f}\n'\
          .format(net, sta, p_arr, s_arr, s_amp, p_snr, s_snr))

