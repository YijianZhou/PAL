"""
params
    data_dir_
    time range
return
    data_dict - 3 chn data paths for each sta
"""
import os, glob
import numpy as np

def get_zsy(data_dir, datetime):
    """ get data paths (in dict) from dir, for certain date
    data paths for ZSY network:
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

def get_sta_dict(sta_file):
    """ get station dict, given sta file name (str)
    """
    sta_dict = []
    f = open(sta_file); lines = f.readlines(); f.close()
    for line in lines:
        sta, lon, lat, ele = line.split('\t')
        sta_dict.append((sta, float(lon), float(lat), float(ele)))
    # convert to struct np.array
    sta_dict = np.array(sta_dict, dtype=[('station','O'),
                                         ('longitude','O'),
                                         ('latitude','O'),
                                         ('elevation','O')])
    return sta_dict
