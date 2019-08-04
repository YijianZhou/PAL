fname = '/home/zhouyj/California/preprocess/station.dat'
fout  = 'input/rc.sta'
f=open(fname); lines=f.readlines(); f.close()
out=open(fout,'w')
for line in lines:
    net, sta, chn, lon, lat, ele = line.split(',')
    lon = abs(float(lon))
    lat = abs(float(lat))
    ele = int(ele)
    lat_deg = int(lat)
    lat_min = 60*(lat-int(lat))
    lon_deg = int(lon)
    lon_min = 60*(lon-int(lon))
    out.write("{:<5} {}  {}  {} {:7.4f}N{} {:7.4f}W{:4}0.2     0.00  0.00  0.00  0.00 3  0.00--HHZ \n"\
        .format(sta, net[-2:], chn, lat_deg, lat_min, lon_deg, lon_min, ele))
out.close()
