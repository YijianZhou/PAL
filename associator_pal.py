import numpy as np

class PSP_Assoc(object):
  """ Associate P- & S-pick Pairs (PSP) by searching ot and loc clustering
  Inputs
    sta_dict: station location dict
    xy_margin: ratio of lateral (x-y) margin relative to the station range
    xy_grid: grid width for x-y axis (in degree)
    z_grids: grids for z axis (in km)
    ot_dev: max time dev for ot assoc
    max_res: threshold for P travel time res
    max_drop: each pick can only be dropped for max_drop times before being associated 
    min_sta: min number of station to alert a detection
    *note: lateral distance (x-y) in degree; depth in km; elevation in m
  Usage
    import associator_pal
    associator = associator_pal.PSP_Assoc(sta_dict)
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
               max_drop  = 1, 
               min_sta   = 4):
    self.sta_dict = sta_dict
    self.xy_margin = xy_margin
    self.xy_grid = xy_grid
    self.z_grids = z_grids
    self.vp = vp
    self.ot_dev = ot_dev
    self.max_res = max_res
    self.max_drop = max_drop
    self.min_sta = min_sta
    self.tt_dict = self.calc_tt()

  def associate(self, picks, out_ctlg=None, out_pha=None):
    events_loc, events_pick = [], []
    num_picks = len(picks)
    if num_picks==0: return 
    picks = np.sort(picks, order='sta_ot')
    # calc num of ot neighbors 
    num_nbr = np.zeros(num_picks) 
    num_drop = np.zeros(num_picks) 
    for ii in range(num_picks): 
        num_nbr[ii] = sum(abs(picks['sta_ot']-picks['sta_ot'][ii]) < self.ot_dev)
    # assoc each cluster
    print('-'*40+'\n'+'detected events:')
    for _ in range(num_picks):
        if np.amax(num_nbr) < self.min_sta: break 
        # 1. ot assoc
        ots = picks['sta_ot']
        ot_i = ots[np.argmax(num_nbr)]
        to_assoc_idx = np.where(abs(ots-ot_i) < self.ot_dev)[0]
        # 2. loc assoc
        event_loc, event_pick, assoc_idx, drop_idx = self.assoc_loc(picks[to_assoc_idx])
        if len(event_loc)>0: 
            # 3. calc mag
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
            events_loc.append(event_loc_mag)
            event_pick.append(event_pick)
        # del picks that are associated or dropped to many times
        drop_idx = np.array(drop_idx, dtype=np.int32) + to_assoc_idx[0]
        assoc_idx = np.array(assoc_idx, dtype=np.int32) + to_assoc_idx[0]
        num_drop[drop_idx] += 1
        to_del = np.unique(np.concatenate([assoc_idx, np.where(num_drop > self.max_drop)[0]]))
        # renew num_nbr 
        to_renew_idx = np.where(abs(ots-ot_i) < 2*self.ot_dev)[0]
        for idx in to_renew_idx:
            num_nbr[idx] -= sum(abs(ots[to_del]-ots[idx]) < self.ot_dev)
        # update picks, num_nbr, and num_drop
        picks = np.delete(picks, to_del)
        num_nbr = np.delete(num_nbr, to_del)
        num_drop = np.delete(num_drop, to_del)
    if not out_ctlg or not out_pha: return events_loc, events_pick 
    else: return

  def assoc_loc(self, picks):
    res_ttp_mat = 0 # P travel time res
    num_sta_mat = 0 # number of associated stations
    det_dict = {} # potential det loc for each sta
    ot = picks['sta_ot'][len(picks)//2]
    for ii, pick in enumerate(picks):
        net_sta = pick['net_sta']
        pick_key = '%s_%s'%(ii, net_sta)
        ttp_obs = pick['tp'] - ot # pick time to travel time
        ttp_pred = self.tt_dict[net_sta]
        res_i = abs(ttp_pred - ttp_obs)
        is_det = (res_i < self.max_res).astype(float)
        res_i[res_i >= self.max_res] = 0.
        # update res_mat, num_sta_mat, and det_dict
        res_ttp_mat += res_i
        num_sta_mat += is_det
        det_dict[pick_key] = is_det
    # find loc of min res (grid search location)
    num_sta = np.amax(num_sta_mat)
    if num_sta < self.min_sta: return [],[],[],list(range(len(picks)))
    res_ttp_mat /= num_sta
    res_ttp_mat [num_sta_mat < num_sta] = np.inf
    res = np.amin(res_ttp_mat)
    zi, xi, yi = np.unravel_index(np.argmin(res_ttp_mat), res_ttp_mat.shape)
    lon = self.lon_min + xi * self.xy_grid
    lat = self.lat_min + yi * self.xy_grid
    dep = self.z_grids[zi]
    # find associated phase & index
    event_pick, assoc_idx, drop_idx = [], [], []
    for ii,pick in enumerate(picks):
        if det_dict['%s_%s'%(ii,pick['net_sta'])][zi,xi,yi]==1.:
            event_pick.append(pick)
            assoc_idx.append(ii)
        else: drop_idx.append(ii)
    # output as dict
    event_loc = {'evt_ot' : ot, 
                 'evt_lon': round(lon,2), 
                 'evt_lat': round(lat,2),
                 'evt_dep': round(dep,0),
                 'res': round(res,1)}
    return event_loc, event_pick, assoc_idx, drop_idx

  # calc P travel time table
  def calc_tt(self):
    print('making time table')
    tt_dict = {}
    # get x-y range: sta range + margin
    sta_loc = self.sta_dict.values()
    lat = [sta_loc[0] for sta_loc in self.sta_dict.values()]
    lon = [sta_loc[1] for sta_loc in self.sta_dict.values()]
    lon_margin = self.xy_margin * (np.amax(lon) - np.amin(lon))
    lat_margin = self.xy_margin * (np.amax(lat) - np.amin(lat))
    lon_min, lon_max = np.amin(lon)-lon_margin, np.amax(lon)+lon_margin
    lat_min, lat_max = np.amin(lat)-lat_margin, np.amax(lat)+lat_margin
    cos_lat = np.cos((lat_max+lat_min)/2 * np.pi/180)
    # set x-y grid
    x_num = int((lon_max-lon_min) / self.xy_grid)
    y_num = int((lat_max-lat_min) / self.xy_grid)
    # calc P travel time table
    for net_sta, [sta_lat,sta_lon,sta_ele,_] in self.sta_dict.items():
        ttp = -np.ones([len(self.z_grids), x_num, y_num])
        for xi in range(x_num):
          for yi in range(y_num):
            for zi,dep in enumerate(self.z_grids):
                dx = 111 * (lon_min + xi*self.xy_grid - sta_lon) * cos_lat
                dy = 111 * (lat_min + yi*self.xy_grid - sta_lat)
                dz = dep + sta_ele/1000.
                ttp[zi,xi,yi] = np.sqrt(dx**2 + dy**2 + dz**2) / self.vp
        tt_dict[net_sta] = ttp
    self.lat_min, self.lon_min = lat_min, lon_min
    return tt_dict

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

  def write_catalog(self, event_loc, out_ctlg):
    ot  = event_loc['evt_ot']
    lon = event_loc['evt_lon']
    lat = event_loc['evt_lat']
    dep = event_loc['evt_dep']
    mag = event_loc['mag'] if 'mag' in event_loc else -1
    out_ctlg.write('{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag))

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

