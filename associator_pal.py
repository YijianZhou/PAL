import numpy as np

class TS_Assoc(object):

  """ Associate picks by searching ot (time, T) and loc (space, S) cluster
  Inputs
    sta_dict: station location dict
    edge_width: ratio of edges relative to the station range
    xy_grid: grid width for x-y axis (in degree)
    ot_dev: time dev for ot assoc
    max_res: threshold for P travel time res
    assoc_num: min num to alert a detection
    *note: lateral spatial range (x-y) in degree; depth in km; elevation in m
  Usage
    import associator_pal
    associator = associator_pal.TS_Assoc(sta_dict)
    associator.associate(picks, out_ctlg, out_pha)
  """
  
  def __init__(self,
               sta_dict,
               edge_width = 0.2,
               xy_grid    = 0.02,
               vp         = 5.9,
               ot_dev     = 3.,
               max_res    = 2.,
               assoc_num  = 4):

    self.sta_dict   = sta_dict
    self.edge_width = edge_width
    self.xy_grid    = xy_grid
    self.vp         = vp
    self.ot_dev     = ot_dev
    self.max_res    = max_res
    self.assoc_num  = assoc_num
    self.tt_dict    = self.calc_tt()


  def associate(self, picks, out_ctlg=None, out_pha=None):
    # 1. temporal assoc: picks --> event_picks
    event_picks = self.assoc_ot(picks)
    # 2. spatial assoc: event_pick --> event_loc
    for event_pick in event_picks:
        event_loc, event_pick = self.assoc_loc(event_pick)
        if len(event_loc)==0: continue
        # 3. estimate magnitude
        event_loc_mag = self.calc_mag(event_pick, event_loc)
        # write catalog and phase
        if out_ctlg: self.write_catalog(event_loc_mag, out_ctlg)
        if out_pha: self.write_phase(event_loc_mag, event_pick, out_pha)
        if not out_ctlg and not out_pha: return event_loc_mag


  # 1. temporal assoc by ot clustering: picks --> event_picks
  def assoc_ot(self, picks):
    event_picks = []
    num_picks = len(picks)
    if num_picks==0: print('no event detected'); return []
    picks = np.sort(picks, order='sta_ot')
    # calc num of ot neighbors (num_nbr)
    num_nbr = np.zeros(num_picks)
    ots = picks['sta_ot']
    for i in range(len(ots)):
        is_nbr = abs(ots-ots[i]) < self.ot_dev
        num_nbr[i] = sum(is_nbr.astype(float))
    # assoc each cluster
    for _ in range(num_picks):
        if np.amax(num_nbr) < self.assoc_num: break
        # ot assoc
        ot_i = ots[np.argmax(num_nbr)]
        to_assoc = abs(ots-ot_i) < self.ot_dev
        event_picks.append(picks[to_assoc])
        num_nbr[to_assoc] = 0
        # renew num_nbr
        to_renew = (abs(ots-ot_i)>self.ot_dev)*(abs(ots-ot_i)<2*self.ot_dev)
        ots_todel = ots[to_assoc]
        nbr_todel = [sum(abs(ots_todel-ot_j)<self.ot_dev) for ot_j in ots[to_renew]]
        num_nbr[to_renew] -= nbr_todel
    print('origin time clustering: %s event candidates'%len(event_picks))
    return event_picks


  # 2. spatial assoc by locate event: event_pick --> event_loc
  def assoc_loc(self, event_pick):
    res_ttp_mat = 0 # P travel time res
    num_sta_mat = 0 # number of associated stations
    det_dict = {} # potential det loc for each sta
    ot = event_pick['sta_ot'][len(event_pick)//2]
    for pick in event_pick:
        net_sta = pick['net_sta']
        ttp_obs = pick['tp'] - ot # pick time to travel time
        ttp_pred = self.tt_dict[net_sta]
        res_i = abs(ttp_pred - ttp_obs)
        is_det = (res_i <= self.max_res).astype(float)
        res_i[res_i > self.max_res] = 0.
        # update res_mat, det_dict, and num_sta_mat
        if np.amax(is_det)==0: continue
        res_ttp_mat += res_i
        det_dict[net_sta] = is_det
        num_sta_mat += is_det
    # find loc of min res (grid search location)
    num_sta = np.amax(num_sta_mat)
    if num_sta < self.assoc_num: return [],[]
    res_ttp_mat /= num_sta
    res_ttp_mat[num_sta_mat != num_sta] = np.inf
    res = np.amin(res_ttp_mat)
    x, y = np.unravel_index(np.argmin(res_ttp_mat), res_ttp_mat.shape)
    lon = self.lon_range[0] + x * self.x_grid
    lat = self.lat_range[0] + y * self.y_grid
    # find associated phase
    event_pick = [pick for pick in event_pick if pick['net_sta'] in det_dict and det_dict[pick['net_sta']][x][y] == 1.]
    # output
    print('locate event: {} {:.2f} {:.2f} | res {:.1f}s'.format(ot, lat, lon, res))
    event_loc = {'evt_ot' : ot, 
                 'evt_lon': round(lon,2), 
                 'evt_lat': round(lat,2),
                 'res': round(res,1)}
    return event_loc, event_pick


  # calc time table
  def calc_tt(self):
    print('making time table')
    tt_dict = {}
    dist_grid = 111 * self.xy_grid # in km
    # get x-y range: sta range + edge width
    sta_loc = np.array(list(self.sta_dict.values()))
    lat, lon = sta_loc['sta_lat'], sta_loc['sta_lon']
    # calc edge width
    lon_edge = self.edge_width * (np.amax(lon) - np.amin(lon))
    lat_edge = self.edge_width * (np.amax(lat) - np.amin(lat))
    # get range
    lon_range = [np.amin(lon) - lon_edge, np.amax(lon) + lon_edge]
    lat_range = [np.amin(lat) - lat_edge, np.amax(lat) + lat_edge]
    # set x-y axis
    cos_lat = np.cos(np.mean(lat_range) * np.pi/180)
    x_range = range(int((lon_range[1]-lon_range[0])*cos_lat / self.xy_grid))
    y_range = range(int((lat_range[1]-lat_range[0]) / self.xy_grid))
    # calc time table
    for net_sta, sta_loc in self.sta_dict.items():
        # convert to x-y axis
        lon = sta_loc['sta_lon']
        lat = sta_loc['sta_lat']
        ele = sta_loc['sta_ele'] / 1000. # m to km
        sta_x = int((lon-lon_range[0])*cos_lat / self.xy_grid)
        sta_y = int((lat-lat_range[0]) / self.xy_grid)
        # calc P travel time
        ttp = -np.ones([len(x_range), len(y_range)])
        for i,x in enumerate(x_range):
          for j,y in enumerate(y_range):
            dist = dist_grid * np.sqrt((x-sta_x)**2 + (y-sta_y)**2)
            dep  = ele + 5 # init event depth 5km
            ttp[i,j] = np.sqrt(dep**2 + dist**2) / self.vp
        tt_dict[net_sta] = ttp
    self.lon_range = lon_range
    self.lat_range = lat_range
    self.x_grid = self.xy_grid * cos_lat
    self.y_grid = self.xy_grid
    return tt_dict


  # calc mag with picks (s_amp)
  def calc_mag(self, event_pick, event_loc, event_dep=5):
    num_sta = len(event_pick)
    mag = -np.ones(num_sta)
    for i,pick in enumerate(event_pick):
        # get sta_loc
        sta_loc = self.sta_dict[pick['net_sta']]
        # get S amp
        if 's_amp' not in pick.dtype.names: continue
        amp = pick['s_amp'] * 1e6 # m to miu m
        # calc epi dist
        dist_lat = 111*(sta_loc['sta_lat'] - event_loc['evt_lat'])
        dist_lon = 111*(sta_loc['sta_lon'] - event_loc['evt_lon']) \
                   * np.cos(sta_loc['sta_lat'] * np.pi/180)
        dist = np.sqrt(dist_lon**2 + dist_lat**2 + event_dep**2)
        mag[i] = np.log10(amp) + np.log10(dist)
    # remove one outlier
    mag_dev = abs(mag - np.median(mag))
    mag = np.delete(mag, np.argmax(mag_dev))
    event_loc['mag'] = round(np.median(mag),2)
    print('estimated magnitude: {:.1f} +- {:.1f}'.format(np.median(mag), np.std(mag)))
    return event_loc


  # write event loc into catalog
  def write_catalog(self, event_loc, out_ctlg):
    ot  = event_loc['evt_ot']
    lon = event_loc['evt_lon']
    lat = event_loc['evt_lat']
    mag = event_loc['mag'] if 'mag' in event_loc else -1
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
        net_sta = pick['net_sta']
        tp = pick['tp']
        ts = pick['ts']
        s_amp = pick['s_amp'] if 's_amp' in pick.dtype.names else -1
        p_snr = pick['p_snr'] if 'p_snr' in pick.dtype.names else -1
        out_pha.write('{},{},{},{},{:.1f}\n'.format(net_sta, tp, ts, s_amp, p_snr))

