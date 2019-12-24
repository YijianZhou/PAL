""" Data Pipeline: Interface for PAD I/O
"""
import os, glob
import numpy as np
from obspy import read, UTCDateTime
import multiprocessing as mp


def get_data(data_dir, date, apply_mp=False):
    """ get data dict
    Input
        data_dir (str): root dir, e.g. root/net/sta/yyyy/mm/dd/net.sta.yyyymmdd.chn.sac
        date (obspy.UTCDateTime): which day of data to get
    Output
        stream_paths = data_dict[net_sta]
        (note: use net.sta to seperate sta from different net)
    """
    # get data paths
    print('reading data: {}'.format(date.date))
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

    # read data
    if apply_mp:
        pool = mp.Pool(processes=5)
        outputs = pool.starmap(read_st, data_dict.items())
        pool.close()
        pool.join()
        for net_sta, st in outputs: data_dict[net_sta] = st
    else:
        for net_sta, st_paths in data_dict.items(): 
            data_dict[net_sta] = read_st(net_sta, st_paths)[1]
    return data_dict


def read_st(net_sta, st_paths):
    # read data
    st  = read(st_paths[0])
    st += read(st_paths[1])
    st += read(st_paths[2])
    # change header
    for i in range(3): 
        st[i].stats.network, st[i].stats.station = net_sta.split('.')
    return net_sta, st


def get_sta(sta_file):
    """ get station dict
    Input
        sta_file: path for sta file, e.g. net, sta, stla, stlo, stel
        (note: use net.sta to seperate sta from different net)
    Output: sta_dict
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


def get_picks(ppk_dir, date):
    """ get picks
    """
    picks = []
    # set output format
    dtype = [('net','O'),
             ('sta','O'),
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
        net, sta, ot0, tp, ts, amp, p_snr, s_snr, fd = line.split(',')
        ot0 = UTCDateTime(ot0)
        tp  = UTCDateTime(tp)
        ts  = UTCDateTime(ts)
        amp = float(amp)
        p_snr = float(p_snr)
        s_snr = float(s_snr)
        fd = float(fd)
        picks.append((net, sta, ot0, tp, ts, amp, p_snr, s_snr, fd))
    return np.array(picks, dtype=dtype)


""" customized data_pipelines
"""
def get_xls(data_dir, datetime):

    data_dict = {}
    year  = str(datetime.year)
    month = str(datetime.month).zfill(2)
    day   = str(datetime.day).zfill(2)
    data_paths = os.path.join(data_dir, year, month, day, '*')
    data_paths = sorted(glob.glob(data_paths))
    for data_path in data_paths:
        file_name = os.path.split(data_path)[-1]
        net_sta = '.'.join(file_name.split('.')[0:2])
        if net_sta in data_dict: data_dict[net_sta].append(data_path)
        else: data_dict[net_sta] = [data_path]
    # drop bad sta
    todel = [net_sta for net_sta in data_dict if len(data_dict[net_sta])!=3]
    for net_sta in todel: data_dict.pop(net_sta)
    return data_dict


def get_rc(data_dir, datetime):
    """ get data paths (in dict) from dir, for certain date
    data paths for CI network:
        net/sta/year/month/day/[net].[sta].[year].[jday].[chn].SAC
    """
    data_dict = {}
    year  = str(datetime.year)
    month = str(datetime.month).zfill(2)
    day   = str(datetime.day).zfill(2)
    data_paths = os.path.join(data_dir, year, month, day, '*')
    data_paths = sorted(glob.glob(data_paths))
    for data_path in data_paths:
        file_name = os.path.split(data_path)[-1]
        sta = file_name.split('.')[2]
        if sta in data_dict: data_dict[sta].append(data_path)
        else: data_dict[sta] = [data_path]
    # drop bad sta
    todel = [sta for sta in data_dict if len(data_dict[sta])!=3]
    for sta in todel: data_dict.pop(sta)
    return data_dict


def get_rc_sta(sta_file):
    """ get station dict, given sta file name (str)
    """
    sta_dict = {}
    dtype = [('sta_lon','O'), ('sta_lat','O'), ('sta_ele','O')]
    f = open(sta_file); lines = f.readlines(); f.close()
    for line in lines:
        net, sta, chn, lon, lat, ele = line.split(',')
        if net not in sta_dict: sta_dict[net] = {}
        sta_dict[net][sta] = np.array((float(lon),float(lat),float(ele)), dtype=dtype)
    return sta_dict

