"""
params
    data_dir
    datetime: obspy.UTCDateTime
return
    data_dict - 3 chn data paths for each sta
"""
import os, glob
import numpy as np
from obspy import UTCDateTime

def get_xj(data_dir, datetime):
    """ get data paths (in dict) from dir, for certain date
    data paths for XJ network:
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
        sta = file_name.split('.')[1]
        if sta in data_dict: data_dict[sta].append(data_path)
        else: data_dict[sta] = [data_path]
    # drop bad sta
    todel = [sta for sta in data_dict if len(data_dict[sta])!=3]
    for sta in todel: data_dict.pop(sta)
    return data_dict


def get_ci(data_dir, datetime):
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


def get_xj_picks(ppk_dir, datetime):
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
             ('freq_domn','O')]
    year  = str(datetime.year)
    month = str(datetime.month).zfill(2)
    day   = str(datetime.day).zfill(2)
    fname = '-'.join([year, month, day]) + '.ppk'
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


def get_sta_dict(sta_file):
    """ get station dict, given sta file name (str)
    """
    sta_dict = {}
    dtype = [('sta_lon','O'), ('sta_lat','O'), ('sta_ele','O')]
    f = open(sta_file); lines = f.readlines(); f.close()
    for line in lines:
        net, sta, lon, lat, ele = line.split('\t')
        if net not in sta_dict: sta_dict[net] = {}
        sta_dict[net][sta] = np.array((float(lon),float(lat),float(ele)), dtype=dtype)
    return sta_dict


def get_ci_sta(sta_file):
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
