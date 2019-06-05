import numpy as np

class TS_Det(object):

  """ Detect by ot (time, T) and ttp (space, S) cluater
  Input
    sta_dict: station location info
    resp_dict: instrumental response dict for all networks (cnt/(m/s))
    side_width: ratio of sides relative to the station range
    xy_grid: grid width for x-y axis (in degree)
    vp, vs: P&S velo of the uniform model
    ot_dev: time dev for ot assoc
    ttp_dev: time dev for ttp assoc
    assoc_num: min num to alert a detection
    *lateral spatial range (x-y) in degree; depth in km; elevation in m
  Usage
    import detectors
    detector = detectors.Simple_Det(sta_dict, resp_dict)
    event_picks = detector.pick2event(picks)
    for event_pick in event_picks:
        event_loc, event_pick = detector.locate(event_pick)
        if len(event_loc)==0: continue
        event_loc_mag = detector.calc_mag(event_pick, event_loc)
        detector.write_catalog(event_loc_mag, out_ctlg)
        detector.write_phase(event_loc_mag, event_pick, out_pha)
  """
  
  def __init__(self,
               sta_dict,
               resp_dict,
               side_width = 0.2,
               xy_grid    = 0.02,
               vp         = 5.9,
               vs         = 3.4,
               ot_dev     = 3.,
               ttp_dev    = 2.,
               assoc_num  = 4):

    self.sta_dict   = sta_dict
    self.resp_dict  = resp_dict
    self.side_width = side_width
    self.xy_grid    = xy_grid
    self.vp         = vp
    self.vs         = vs
    self.ot_dev     = ot_dev
    self.ttp_dev    = ttp_dev
    self.assoc_num  = assoc_num
    self.time_table = self.calc_tt(sta_dict)


  def pick2event(self, picks):
    event_picks = []
    num_picks = len(picks)
    picks = np.sort(picks, order='org_t0')

    # assoc each cluster
    for _ in range(num_picks):
        # calc neighbor num
        ots = picks['org_t0']
        if len(ots)==0: break
        num_nbr = np.zeros(len(ots))
        for i, oti in enumerate(ots):
            is_nbr = abs(ots-oti) < self.ot_dev
            num_nbr[i] = sum(is_nbr.astype(float))
        if np.amax(num_nbr) < self.assoc_num: break
        oti = ots[np.argmax(num_nbr)]
        is_nbr = abs(ots-oti) < self.ot_dev
        event_pick = picks[is_nbr]
        picks = picks[~is_nbr]
        event_picks.append(event_pick)

    print('associated {} events'.format(len(event_picks)))
    return event_picks


  def locate(self, event_pick):
    """ locate events by grid search
    Input
        event_pick: picks for an event
    Output
        location: lon, lat, dep
        quality: res, num_sta
    Algorithm:
        grid search for a 3d time table
    """
    res_ttp = 0 # P travel time res
    tp_num = 0 # true positive
    tp_dict = {}
    num_sta = len(event_pick)
    ot = event_pick['org_t0'][num_sta//2]
    for pick in event_pick:
        sta = pick['station']
        ttp_obs = pick['p_arr'] - ot # to travel time
        ttp_pred = self.time_table[sta]
        resi = abs(ttp_pred - ttp_obs)
        resi[resi > self.ttp_dev] = 0.
        res_ttp += resi
        is_tp = (resi>0.).astype(float)
        tp_dict[sta] = is_tp
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
    event_pick = [pick for pick in event_pick\
                  if tp_dict[pick['station']][x][y] == 1.]

    # output
    print('locate event: org_time {}, lon {:.2f}, lat {:.2f}, res {:.1f}'\
        .format(ot, lon, lat, min_res))
    event_loc = {'org_time' : ot, 
                 'longitude': round(lon,2), 
                 'latitude': round(lat,2),
                 'residual': round(min_res,1)}
    return event_loc, event_pick


  def calc_tt(self, sta_dict):
    """ calc Time Table for all stations
    Calc 3D travel time table
    x-y-z axis: x - lon, y - lat, z - depth
    Note that x-y range for loc are determined in this step
    """
    # set up
    print('Making time table')
    time_table = {}
    dist_grid = 111 * self.xy_grid # in km

    # get x-y range: sta range + side width
    lon = sta_dict['longitude']
    lat = sta_dict['latitude']
    # calc lon side
    lon_len  = np.amax(lon) - np.amin(lon)
    lon_side = self.side_width * lon_len
    # calc lat side
    lat_len  = np.amax(lat) - np.amin(lat)
    lat_side = self.side_width * lat_len
    # get range
    lon_rng = [np.amin(lon) - lon_side, np.amax(lon) + lon_side]
    lat_rng = [np.amin(lat) - lat_side, np.amax(lat) + lat_side]
    # set x-y axis
    x_rng = range(int((lon_rng[1] - lon_rng[0]) / self.xy_grid))
    y_rng = range(int((lat_rng[1] - lat_rng[0]) / self.xy_grid))

    # calc time table
    for sta_info in sta_dict:
        # convert to x-y axis
        sta = sta_info['station']
        lon = sta_info['longitude']
        lat = sta_info['latitude']
        ele = sta_info['elevation'] / 1000. # m to km
        sta_x = int((lon - lon_rng[0]) / self.xy_grid)
        sta_y = int((lat - lat_rng[0]) / self.xy_grid)

        # calc P travel time
        tt_p = -np.ones([len(x_rng), len(y_rng)])
        for i,x in enumerate(x_rng):
          for j,y in enumerate(y_rng):
            dist = dist_grid * np.sqrt((x-sta_x)**2 + (y-sta_y)**2)
            dep  = ele + 5 # init event depth 5km
            tt_p[i,j] = np.sqrt(dep**2 + dist**2) / self.vp
        time_table[sta] = tt_p

    self.lon_rng = lon_rng
    self.lat_rng = lat_rng
    return time_table


  def calc_mag(self, event_pick, event_loc):
    """ calc magnitude with event picks (amp) and loc
    """
    num_sta = len(event_pick)
    mag  = np.zeros(num_sta)
    dist = np.zeros(num_sta)
    for i,pick in enumerate(event_pick):
        # get S amplitude
        resp = self.resp_dict[pick['network']]
        amp = pick['s_amp'] /resp *1e6 # in miu m
        # get epicentral distance
        dist_lon = 111*(pick['sta_lon'] - event_loc['longitude'])
        dist_lat = 111*(pick['sta_lat'] - event_loc['latitude'])
        dist[i] = np.sqrt(dist_lon**2 + dist_lat**2) # in km
        mag[i]  = np.log10(amp) + np.log10(dist[i])
    mag = mag[dist<2*np.amin(dist)] # drop bad data
    event_loc['magnitude'] = round(np.mean(mag),1)
    print('estimated mag: {:.1f} delta {:.1f}'.\
          format(np.mean(mag), np.std(mag)))
    return event_loc


  def write_catalog(self, event_loc, out_ctlg):
    """ write event loc into catalog
    """
    ot  = event_loc['org_time']
    lon = event_loc['longitude']
    lat = event_loc['latitude']
    mag = event_loc['magnitude']
    res = event_loc['residual']
    out_ctlg.write('{},{},{},5,{},{}\n'\
           .format(ot, lon, lat, mag, res))


  def write_phase(self, event_loc, event_pick, out_pha):
    """ write sta phase into phase file
    """
    ot  = event_loc['org_time']
    lon = event_loc['longitude']
    lat = event_loc['latitude']
    mag = event_loc['magnitude']
    res = event_loc['residual']
    out_pha.write('{},{},{},5,{},{}\n'\
          .format(ot, lon, lat, mag, res))
    for pick in event_pick:
        net   = pick['network']
        sta   = pick['station']
        p_arr = pick['p_arr']
        s_arr = pick['s_arr']
        s_amp = pick['s_amp']
        p_snr = pick['p_snr']
        s_snr = pick['s_snr']
        out_pha.write('{},{},{},{},{:.1f},{:.1f},{:.1f}\n'\
                    .format(net, sta, p_arr, s_arr, s_amp, p_snr, s_snr))

