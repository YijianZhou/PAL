import numpy as np

class Simple_Loc(object):

  """Simple location by grid search
  Input
    sta_dict: station location info
    side_width: ratio of sides relative to the station range
    dep_rng: range of location depth (in km)
    xy_grid: grid width for x-y axis (in degree)
    z_grid: grid width for z axis (in km)
    vp, vs: P&S velo of the uniform model
    *lateral spatial range (x-y) in degree; depth in km; elevation in m
  Usage
    import locators
    locator = locators.Simple_Loc(sta_dict)
    loc = locator.locate(event_pick)
  """
  
  def __init__(self,
               sta_dict,
               side_width = 0.2,
               dep_rng    = [0.,30.],
               xy_grid    = 0.05,
               z_grid     = 2.5,
               vp         = 5.8,
               vs         = 3.2):
               
    self.sta_dict   = sta_dict
    self.side_width = side_width
    self.dep_rng    = dep_rng
    self.xy_grid    = xy_grid
    self.z_grid     = z_grid
    self.vp         = vp
    self.vs         = vs
    self.time_table = self.calc_tt(sta_dict)


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
    res= 0
    num_sta = len(event_pick)
    ot = event_pick['org_time0'][int(num_sta/2)]
    for pick in event_pick:
        sta = pick['station']
        tp = pick['p_arr'] - ot # to travel time
        ts = pick['s_arr'] - ot
        tt_p = self.time_table[sta][0]
        tt_s = self.time_table[sta][1]
        res_p = tt_p - tp
        res_s = tt_s - ts
        res += np.abs(res_p) + np.abs(res_s)

    # find min res in the map
    min_res = np.amin(res)
    x, y, z = np.unravel_index(np.argmin(res), res.shape)
    lon = self.lon_rng[0] + x*self.xy_grid
    lat = self.lat_rng[1] + y*self.xy_grid
    dep = self.dep_rng[2] + z*self.z_grid

    # output
    print('locate event: ot {}, lon {}, lat {}, dep {}km, res ={}'\
          .format(ot, lon, lat, dep, min_res))
    event_loc = {'org_time' : ot, 
                 'longitude': lon, 
                 'latitude' : lat, 
                 'depth'    : dep}
    return event_loc


  def calc_tt(self):
    """ calc Time Table for all stations
    Calc 3D travel time table
    x-y-z axis: x - lon, y - lat, z - depth
    Note that x-y range for loc are determined in this step
    """
    # set up
    time_table = {}
    grid_width = 111*self.xy_grid # in km

    # get x-y range: sta range + side width
    lon = sta_dict['longitude']
    lat = sta_dict['latitude']
    # calc lon side
    lon_len  = np.amax(lon) - np.amin(lon)
    lon_side = self.side_width *lon_len
    # calc lat side
    lat_len  = np.amax(lat) - np.amin(lat)
    lat_side = self.side_width *lat_len
    # get range
    lon_rng = [np.amin(lon)-lon_side, np.amax(lon)+lon_side]
    lat_rng = [np.amin(lat)-lat_side, np.amax(lat)+lat_side]
    dep_rng = self.dep_rng # in km
    # set x-y-z axis
    x_rng = range(int((lon_rng[1] - lon_rng[0]) /self.xy_grid))
    y_rng = range(int((lat_rng[1] - lat_rng[0]) /self.xy_grid))
    z_rng = range(dep_rng[0], dep_rng[1], self.z_grid)

    # calc time table
    for sta_info in sta_dict:
        # convert to x-y axis
        sta = sta_info['station']
        lon = sta_info['longitude']
        lat = sta_info['latitude']
        ele = sta_info['elevation'] /1000 # m to km
        sta_x = int((lon - lon_rng[0]) /self.xy_grid)
        sta_y = int((lat - lat_rng[0]) /self.xy_grid)

        # calc P and S travel time
        tt_p = -np.ones(len(x_rng), len(y_rng), len(z_rng))
        tt_s = -np.ones(len(x_rng), len(y_rng), len(z_rng))
        for i,x in enumerate(x_rng):
          for j,y in enumerate(y_rng):
            for k,z = in enumerate(z_rng):
                dist = grid_width *np.sqrt((x-sta_x)**2 + (y-sta_y)**2)
                tt_p[i,j,k] = np.sqrt((ele+z)**2 + dist**2) /self.vp
                tt_s[i,j,k] = np.sqrt((ele+z)**2 + dist**2) /self.vs
        time_table[sta] = [tt_p, tt_s]

    self.lon_rng = lon_rng
    self.lat_rng = lat_rng
    return time_table


  def calc_mag(self, event_pick, event_loc):
    """ calc magnitude with event picks (amp) and loc
    """
    mag = []
    for pick in event_pick:
        dist_lon = 111*(pick['sta_lon'] - event_loc['longitude'])
        dist_lat = 111*(pick['sta_lat'] - event_loc['latitude'])
        dist = np.sqrt(dist_lon**2 + dist_lat**2) # in km
        magi = 
        mag.append(magi)
    mag = np.array(mag)
    event_loc['magnitude'] = np.mean(mag)
    print('estimated mag: {} ({})'.format(np.mean(mag), np.std(mag))
    return event_loc


  def write(self, event_loc, out_ctlg):
    """ write event loc into catalog
    format: 
    """
    ot  = event_loc['org_time']
    lon = event_loc['longitude']
    lat = event_loc['latitude']
    dep = event_loc['depth']
    mag = event_loc['magnitude']
    out_ctlg.write('{},{},{},{},{}\n'\
                  .format(ot, lon, lat, dep, mag))
    return

