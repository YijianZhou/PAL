""" Make input station file for HypoDD
"""
import os
import config

# i/o paths
cfg = config.Config()
fout=open('input/station.dat','w')
done_list = []

f=open(cfg.fsta); lines=f.readlines(); f.close()
for line in lines:
    codes = line.split(',')
    net, sta = codes[0].split('.')
    if sta in done_list: continue
    lat = float(codes[1])
    lon = float(codes[2])
    fout.write('{} {} {}\n'.format(sta, lat, lon))
    done_list.append(sta)
fout.close()
