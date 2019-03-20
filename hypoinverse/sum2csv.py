fname = 'xj.sum'
fout  = 'xj.csv'
f=open(fname); lines=f.readlines(); f.close()
out=open(fout,'w')

for line in lines:
    info = line.split()
    date = info[0]
    hrmn = info[1]
    sec  = info[2]
    dtime = date+hrmn+sec.zfill(5)
    lat_deg = float(info[3])
    lat_min = float(info[4])
    lat = lat_deg + lat_min/60
    lon_deg = float(line.split('E')[0].split()[-1])
    lon_min = float(line.split('E')[1].split()[0])
    lon = lon_deg + lon_min/60
    dep = float(line.split('E')[1].split()[1])
    out.write('{},{:.4f},{:.4f},{:.1f}\n'.format(dtime, lat, lon, dep))

