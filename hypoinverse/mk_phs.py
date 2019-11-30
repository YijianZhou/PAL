""" make phase input file for HypoInverse (COP 3 format)
  change net to choose input phase file
"""
from obspy import UTCDateTime

net = 'zsy'
fpha  = '../output/%s/phase_%s.dat'%(net,net)
fout  = 'input/%s.phs'%net
lat_code, lon_code = 'N', 'E'
f=open(fpha); lines =f.readlines(); f.close()
out=open(fout,'w')
mag_corr = 2. # hypoInv do not support neg mag

def split_datetime(dtime):
    yr  = dtime.year
    mon = dtime.month
    day = dtime.day
    date = '{}{:0>2}{:0>2}'.format(yr, mon, day)
    hr  = dtime.hour
    mi  = dtime.minute
    sec = dtime.second
    msc = int(dtime.microsecond/1e4)
    time = '{:0>2}{:0>2}{:0>2}{:0>2}'.format(hr, mi, sec, msc)
    return date, time

idx=0
for line in lines:
  codes = line.split(',')
  if len(codes)==5:
    # write head line
    ot, lat, lon, mag = codes[0:4]
    ot = UTCDateTime(ot)
    date, time = split_datetime(ot)
    mag = max(float(mag) + mag_corr, 0.)
    lat = abs(float(lat))
    lon = abs(float(lon))
    lon_deg = int(lon)
    lon_min = int(100*60*(lon-int(lon)))
    lat_deg = int(lat)
    lat_min = int(100*60*(lat-int(lat)))
    lat = '{:0>2}{}{:0>4}'.format(lat_deg, lat_code, lat_min)
    lon = '{:0>3}{}{:0>4}'.format(lon_deg, lon_code, lon_min)
    if idx!=0: out.write('\n')
    out.write('{}{}{} {}L{:3.2f}{}{:>10}L\n'\
      .format(date+time, lat, lon, ' '*90, mag, ' '*9, idx))
    idx+=1
  else:
    # write sta line
    net, sta, tp, ts  = codes[0:4]
    tp = UTCDateTime(tp)
    ts = UTCDateTime(ts)
    date = split_datetime(tp)[0]
    ts = ts -(tp - tp.second - tp.microsecond/1e6)
    tp = split_datetime(tp)[1]
    ts = int(100*ts)
    out.write('{:<5}{}  HHZ IPU0{} {}{}{:5}ES 1 \n'\
      .format(sta, net[-2:], date+tp[0:4], tp[4:],' '*7, ts))

