""" Data Pipeline: Interface for PAD I/O
"""
import os, glob
import numpy as np
from obspy import read, UTCDateTime

class Data(object):

  def __init__(self, resp_dict):
    self.resp_dict = resp_dict


  def get_data_dict(self, date, data_dir):
    """ get data path dict
    Inputs
      data_dir (str): root dir, e.g. root/net/sta/yyyy/mm/dd/net.sta.yyyymmdd.chn.sac
      date (obspy.UTCDateTime): which day of data to get
    Outputs
      st_paths = data_dict[net_sta]
      *note: use net.sta to seperate sta from different net
    """
    # get data paths
    data_dict = {}
    date_dir = '{:0>4}/{:0>2}/{:0>2}'.format(date.year, date.month, date.day)
    st_paths = sorted(glob.glob(os.path.join(data_dir, date_dir, '*')))
    for st_path in st_paths:
        fname = os.path.split(st_path)[-1]
        net_sta = '.'.join(fname.split('.')[0:2])
        if net_sta in data_dict: data_dict[net_sta].append(st_path)
        else: data_dict[net_sta] = [st_path]
    # drop bad sta
    todel = [net_sta for net_sta in data_dict if len(data_dict[net_sta])!=3]
    for net_sta in todel: data_dict.pop(net_sta)
    return data_dict


  def read_data(self, st_paths):
    # read data
    print('reading stream: {}'.format(st_paths[0]))
    st  = read(st_paths[0])
    st += read(st_paths[1])
    st += read(st_paths[2])
    # change header
    net, sta = os.path.split(st_paths[0])[-1].split('.')[0:2]
    for i in range(3): 
        st[i].stats.network, st[i].stats.station = net, sta
        st[i].data /= self.resp_dict[net]
    return st


  def get_sta_dict(self, sta_file):
    """ get station location dict
    Inputs
      sta_file: path for sta file, e.g. net, sta, stla, stlo, stel
      *note: use net.sta to seperate sta from different net
    Outputs
      stla, stlo, stel = sta_dict[net_sta]
    """
    sta_dict = {}
    dtype = [('sta_lon','O'), ('sta_lat','O'), ('sta_ele','O')]
    f = open(sta_file); lines = f.readlines(); f.close()
    for line in lines:
        net, sta, lon, lat, ele = line.split('\t')
        net_sta = '.'.join([net, sta])
        sta_dict[net_sta] = np.array((float(lon),float(lat),float(ele)), dtype=dtype)
    return sta_dict


  def get_picks(self, date, ppk_dir):
    """ get PAD picks (for assoc)
    """
    picks = []
    # set output format
    dtype = [('net_sta','O'),
             ('sta_ot','O'),
             ('p_arr','O'),
             ('s_arr','O'),
             ('s_amp','O'),
             ('p_snr','O'),
             ('s_snr','O'),
             ('freq_dmnt','O')]
    fname = str(date.date) + '.ppk'
    ppk_path = os.path.join(ppk_dir, fname)
    f=open(ppk_path); lines=f.readlines(); f.close()
    for line in lines:
        codes = line.split(',')
        net_sta = '.'.join(codes[0:2])
        sta_ot, tp, ts = [UTCDateTime(t) for t in codes[2:5]]
        amp, p_snr, s_snr, fd = [float(x) for x in codes[5:9]]
        picks.append((net_sta, sta_ot, tp, ts, amp, p_snr, s_snr, fd))
    return np.array(picks, dtype=dtype)


""" customized data_pipelines
"""
