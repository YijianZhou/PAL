fname = '/data3/XJ_SAC/header/station.dat'
fout  = 'xj.sta'
f=open(fname); lines=f.readlines(); f.close()
out=open(fout,'w')
chn='HHZ'
for line in lines:
    net, sta, lon, lat, ele = line.split('\t')
    lon = float(lon)
    lat = float(lat)
    ele = int(ele)
    lat_int = int(lat)
    lat_flt = 60*(lat-int(lat))
    lon_int = int(lon)
    lon_flt = 60*(lon-int(lon))
    out.write("{:<5} {}  {}  {} {:7.4f}N{} {:7.4f}E{:4}0.2     0.00  0.00  0.00  0.00 3  0.00--HHZ \n"\
        .format(sta, net[-2:], chn, lat_int, lat_flt, lon_int, lon_flt, ele))
out.close()
