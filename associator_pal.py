import numpy as np

class TS_Assoc(object):

  """ Associate picks by searching ot (time, T) and loc (space, S) clustering
  Inputs
    sta_dict: station location dict
    xy_margin: ratio of lateral (x-y) margin relative to the station range
    xy_grid: grid width for x-y axis (in degree)
    z_grids: grids for z axis (in km)
    ot_dev: max time dev for ot assoc
    max_res: threshold for P travel time res
    min_sta: min number of station to alert a detection
    *note: lateral distance (x-y) in degree; depth in km; elevation in m
  Usage
    import associator_pal
    associator = associator_pal.TS_Assoc(sta_dict)
    associator.associate(picks, out_ctlg, out_pha)
  """
  
  def __init__(self,
               sta_dict,
               xy_margin = 0.2,
               xy_grid   = 0.02,
               z_grids   = [5],
               vp        = 5.9,
               ot_dev    = 2.5,
               max_res   = 1.5,
               min_sta   = 4):

    self.sta_dict  = sta_dict
    self.xy_margin = xy_margin
    self.xy_grid   = xy_grid
    self.z_grids   = z_grids
    self.vp        = vp
    self.ot_dev    = ot_dev
    self.max_res   = max_res
    self.min_sta   = min_sta
    self.tt_dict   = self.calc_tt()


  def associate(self, picks, out_ctlg=None, out_pha=None):
    # 1. temporal assoc: picks --> event_picks
    event_picks = self.assoc_ot(picks)
    events, picks = [], []
    # 2. spatial assoc: event_pick --> event_loc
    print('-'*40)
    print('detected events:')
    for event_pick in event_picks:
        event_loc, event_pick = self.assoc_loc(event_pick)
        if len(event_loc)==0: continue
        # 3. estimate magnitude
        event_loc_mag = self.calc_mag(event_pick, event_loc)
        # screen output
        ot  = event_loc_mag['evt_ot']
        lon = event_loc_mag['evt_lon']
        lat = event_loc_mag['evt_lat']
        dep = event_loc_mag['evt_dep']
        mag = event_loc_mag['mag']
        res = event_loc_mag['res']
        print('{} {} {} {:>2} {} | res {}s'.format(ot, lat, lon, dep, mag, res))
        # write catalog and phase
        if out_ctlg: self.write_catalog(event_loc_mag, out_ctlg)
        if out_pha: self.write_phase(event_loc_mag, event_pick, out_pha)
        events.append(event_loc_mag)
        picks.append(event_pick)
    if not out_ctlg or not out_pha: return events, picks


  # 1. temporal assoc by ot clustering: picks --> event_picks
  def assoc_ot(self, picks):
    event_picks = []
    num_picks = len(picks)
    if num_picks==0: return []
    picks = np.sort(picks, order='sta_ot')
    # calc num of ot neighbors (num_nbr)
    num_nbr = np.zeros(num_picks)
    ots = picks['sta_ot']
    for i in range(len(ots)):
        is_nbr = abs(ots-ots[i]) < self.ot_dev
        num_nbr[i] = sum(is_nbr.astype(float))
    # assoc each cluster
    for _ in range(num_picks):
        if np.amax(num_nbr) < self.min_sta: break
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
    return event_picks


  # 2. spatial assoc by locate event: event_pick --> event_loc
  def assoc_loc(self, event_pick):
    res_ttp_mat = 0 # P travel time res
    num_sta_mat = 0 # number of associated stations
    det_dict = {} # potential det loc for each sta
    ot = event_pick['sta_ot'][len(event_pick)//2]
    bad_idx = []
    for i, pick in enumerate(event_pick):
        net_sta = pick['net_sta']
        ttp_obs = pick['tp'] - ot # pick time to travel time
        ttp_pred = self.tt_dict[net_sta]
        res_i = abs(ttp_pred - ttp_obs)
        is_det = (res_i < self.max_res).astype(float)
        res_i[res_i >= self.max_res] = 0.
        # update res_mat, det_dict, and num_sta_mat
        if np.amax(is_det)==0: bad_idx.append(i); continue
        if net_sta in det_dict: bad_idx.append(i); continue
        res_ttp_mat += res_i
        det_dict[net_sta] = is_det
        num_sta_mat += is_det
    event_pick = np.delete(event_pick, bad_idx)
    # find loc of min res (grid search location)
    num_sta = np.amax(num_sta_mat)
    if num_sta < self.min_sta: return [],[]
    res_ttp_mat /= num_sta
    res_ttp_mat [num_sta_mat < num_sta] = np.inf
    res = np.amin(res_ttp_mat)
    zi, x, y = np.unravel_index(np.argmin(res_ttp_mat), res_ttp_mat.shape)
    lon = self.lon_range[0] + x * self.xy_grid
    lat = self.lat_range[0] + y * self.xy_grid
    dep = self.z_grids[zi]
    # find associated phase
    event_pick = [pick for pick in event_pick if det_dict[pick['net_sta']][zi,x,y]==1.]
    # output
    event_loc = {'evt_ot' : ot, 
                 'evt_lon': round(lon,2), 
                 'evt_lat': round(lat,2),
                 'evt_dep': round(dep,0),
                 'res': round(res,1)}
    return event_loc, event_pick


  # calc time table
  def calc_tt(self):
    print('making time table')
    tt_dict = {}
    # get x-y range: sta range + margin
    sta_loc = self.sta_dict.values()
    lat = [sta_loc[0] for sta_loc in self.sta_dict.values()]
    lon = [sta_loc[1] for sta_loc in self.sta_dict.values()]
    lon_margin = self.xy_margin * (np.amax(lon) - np.amin(lon))
    lat_margin = self.xy_margin * (np.amax(lat) - np.amin(lat))
    lon_range = [np.amin(lon)-lon_margin, np.amax(lon)+lon_margin]
    lat_range = [np.amin(lat)-lat_margin, np.amax(lat)+lat_margin]
    cos_lat = np.cos(np.mean(lat_range) * np.pi/180)
    # set x-y grid
    x_num = int((lon_range[1]-lon_range[0]) / self.xy_grid)
    y_num = int((lat_range[1]-lat_range[0]) / self.xy_grid)
    # calc time table
    for net_sta, [lat,lon,ele,_] in self.sta_dict.items():
        # convert to x-y grid_idx
        sta_x = int((lon-lon_range[0]) / self.xy_grid)
        sta_y = int((lat-lat_range[0]) / self.xy_grid)
        # calc P travel time
        ttp = -np.ones([len(self.z_grids), x_num, y_num])
        for x in range(x_num):
          for y in range(y_num):
            for zi,z in enumerate(self.z_grids):
                dx = 111 * (x-sta_x) * self.xy_grid * cos_lat # degree to km
                dy = 111 * (y-sta_y) * self.xy_grid
                dz = z + ele/1000.
                ttp[zi,x,y] = np.sqrt(dx**2 + dy**2 + dz**2) / self.vp
        tt_dict[net_sta] = ttp
    self.lon_range = lon_range
    self.lat_range = lat_range
    return tt_dict


  # calc mag with picks (s_amp)
  def calc_mag(self, event_pick, event_loc):
    num_sta = len(event_pick)
    mag = -np.ones(num_sta)
    for i,pick in enumerate(event_pick):
        sta_lat, sta_lon, sta_ele = self.sta_dict[pick['net_sta']][0:3]
        # get S amp
        if 's_amp' not in pick.dtype.names: continue
        amp = pick['s_amp'] * 1e6 # m to miu m
        # calc epi dist
        dist_lat = 111*(sta_lat - event_loc['evt_lat'])
        dist_lon = 111*(sta_lon - event_loc['evt_lon']) * np.cos(sta_lat*np.pi/180)
        dist_dep = event_loc['evt_dep'] + sta_ele/1e3
        dist = np.sqrt(dist_lon**2 + dist_lat**2 + dist_dep**2)
        mag[i] = np.log10(amp) + np.log10(dist) + 1
    # remove one outlier
    mag_dev = abs(mag - np.median(mag))
    mag = np.delete(mag, np.argmax(mag_dev))
    event_loc['mag'] = round(np.median(mag),2)
    return event_loc


  # write event loc into catalog
  def write_catalog(self, event_loc, out_ctlg):
    ot  = event_loc['evt_ot']
    lon = event_loc['evt_lon']
    lat = event_loc['evt_lat']
    dep = event_loc['evt_dep']
    mag = event_loc['mag'] if 'mag' in event_loc else -1
    out_ctlg.write('{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag))


  # write sta phase into phase file
  def write_phase(self, event_loc, event_pick, out_pha):
    ot  = event_loc['evt_ot']
    lon = event_loc['evt_lon']
    lat = event_loc['evt_lat']
    dep = event_loc['evt_dep']
    mag = event_loc['mag']
    out_pha.write('{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag))
    for pick in event_pick:
        net_sta = pick['net_sta']
        tp = pick['tp']
        ts = pick['ts']
        s_amp = pick['s_amp'] if 's_amp' in pick.dtype.names else -1
        out_pha.write('{},{},{},{}\n'.format(net_sta, tp, ts, s_amp))

