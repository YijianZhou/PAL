""" Make station input file for hypoInverse
"""
import config

# i/o paths
cfg = config.Config()
fsta = cfg.fsta_in 
fout = cfg.fsta_out
lat_code = cfg.lat_code
lon_code = cfg.lon_code
f=open(fsta); lines=f.readlines(); f.close()
out = open(fout,'w')

for line in lines:
    net, sta, lon, lat, ele = line.split('\t')
    lon = abs(float(lon))
    lat = abs(float(lat))
    ele = int(ele)
    lat_deg = int(lat)
    lat_min = 60*(lat-int(lat))
    lon_deg = int(lon)
    lon_min = 60*(lon-int(lon))
    lat = '{} {:7.4f}{}'.format(lat_deg, lat_min, lat_code)
    lon = '{} {:7.4f}{}'.format(lon_deg, lon_min, lon_code)
    # hypoinverse format 2
    out.write("{:<5} {}  HHZ  {}{}{:4}\n".format(sta, net[-2:], lat, lon, ele))
out.close()
