""" Data Pipeline: interface for data i/o
"""
import os, glob
import numpy as np
import obspy
from obspy import read, UTCDateTime


# get data path dict
def get_data_dict(date, data_dir):
    # get data paths
    data_dict = {}
    net_sta_code = '*/*'
    date_code = '{:0>4}/{:0>2}/{:0>2}'.format(date.year, date.month, date.day)
    st_paths = sorted(glob.glob(os.path.join(data_dir, net_sta_code, date_code, '*')))
    for st_path in st_paths:
        fname = os.path.basename(st_path)
        net_sta = '.'.join(fname.split('.')[0:2])
        if net_sta in data_dict: data_dict[net_sta].append(st_path)
        else: data_dict[net_sta] = [st_path]
    # drop bad sta
    todel = [net_sta for net_sta in data_dict if len(data_dict[net_sta])!=3]
    for net_sta in todel: data_dict.pop(net_sta)
    return data_dict


# read stream data
def read_data(st_paths, sta_dict):
    # read data
    print('reading stream: {}'.format(st_paths[0]))
    try:
        st  = read(st_paths[0])
        st += read(st_paths[1])
        st += read(st_paths[2])
    except: 
        print('bad data!'); return []
    # change header
    net, sta = os.path.basename(st_paths[0]).split('.')[0:2]
    net_sta = '%s.%s'%(net,sta)
    for i in range(3): 
        st[i].stats.network, st[i].stats.station = net, sta
        st[i].data = st[i].data / float(sta_dict[net_sta]['gain'])
    return st


# get station loc & gain dict
def get_sta_dict(sta_file):
    sta_dict = {}
    dtype = [('sta_lat','O'),('sta_lon','O'),('sta_ele','O'),('gain','O')]
    f = open(sta_file); lines = f.readlines(); f.close()
    for line in lines:
        net_sta, lat, lon, ele, gain = line.split(',')
        sta_dict[net_sta] = np.array((float(lat),float(lon),float(ele),float(gain)), dtype=dtype)
    return sta_dict


# get PAL picks (for assoc)
def get_picks(date, pick_dir):
    picks = []
    dtype = [('net_sta','O'),
             ('sta_ot','O'),
             ('tp','O'),
             ('ts','O'),
             ('s_amp','O'),
             ('p_snr','O'),
             ('freq_dom','O')]
    fname = str(date.date) + '.pick'
    pick_path = os.path.join(pick_dir, fname)
    f=open(pick_path); lines=f.readlines(); f.close()
    for line in lines:
        codes = line.split(',')
        net_sta = codes[0]
        sta_ot, tp, ts = [UTCDateTime(t) for t in codes[1:4]]
        amp, p_snr, fd = [float(x) for x in codes[4:7]]
        picks.append((net_sta, sta_ot, tp, ts, amp, p_snr, fd))
    return np.array(picks, dtype=dtype)


""" customized data_pipelines
"""
