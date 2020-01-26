""" Make station input file for HypoDD
"""
import config

# i/o paths
cfg = config.Config()
fsta = cfg.fsta_in
fout = cfg.fsta_out
f=open(fsta); lines=f.readlines(); f.close()
out=open(fout,'w')

for line in lines:
    net, sta, lon, lat, ele = line.split('\t')
    lon = float(lon)
    lat = float(lat)
    out.write('{} {} {}\n'.format(sta, lat, lon))
out.close()
