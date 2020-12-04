""" Make input station file for HypoDD
"""
import os
import config

# i/o paths
cfg = config.Config()
fsta = cfg.fsta_in
fout = cfg.fsta_out
f=open(fsta); lines=f.readlines(); f.close()
fout=open(fout,'w')

for line in lines:
    net_sta, lat, lon, ele, _ = line.split(',')
    net, sta = net_sta.split('.')
    lon = float(lon)
    lat = float(lat)
    fout.write('{} {} {}\n'.format(sta, lat, lon))
fout.close()
