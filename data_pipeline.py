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
    return data_dict


def get_xj_picks(ppk_dir, datetime):
    """ get picks
    """
    picks = []
    # set output format
    dtype = [('network','O'),
             ('station','O'),
             ('sta_lon','O'),
             ('sta_lat','O'),
             ('org_t0','O'),
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
        net, sta, sta_lon, sta_lat, ot0, tp, ts, amp, p_snr, s_snr, fd = line.split(',')
        sta_lon = float(sta_lon)
        sta_lat = float(sta_lat)
        ot0 = UTCDateTime(ot0)
        tp  = UTCDateTime(tp)
        ts  = UTCDateTime(ts)
        amp = float(amp)
        p_snr = float(p_snr)
        s_snr = float(s_snr)
        fd = float(fd)
        picks.append((net, sta, sta_lon, sta_lat, ot0, tp, ts, amp, p_snr, s_snr, fd))
    return np.array(picks, dtype=dtype)


def get_sta_dict(sta_file):
    """ get station dict, given sta file name (str)
    """
    sta_dict = []
    f = open(sta_file); lines = f.readlines(); f.close()
    for line in lines:
        net, sta, lon, lat, ele = line.split('\t')
        sta_dict.append((net, sta, float(lon), float(lat), float(ele)))
    # convert to struct np.array
    sta_dict = np.array(sta_dict, dtype=[('network','O'),
                                         ('station','O'),
                                         ('longitude','O'),
                                         ('latitude','O'),
                                         ('elevation','O')])
    return sta_dict
