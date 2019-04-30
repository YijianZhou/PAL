from obspy import UTCDateTime

net = 'xls'
fpha  = '../output/%s/phase_xj_%s.dat'%(net,net)
fctlg = '../output/%s/catalog_xj_%s.dat'%(net,net)
fout  = 'xj.phs'
f=open(fpha); phas =f.readlines(); f.close()
f=open(fctlg);ctlgs=f.readlines(); f.close()
out=open(fout,'w')

def split_datetime(dtime):
    yr  = str(dtime.year)
    mon = str(dtime.month).zfill(2)
    day = str(dtime.day).zfill(2)
    date = yr + mon + day
    hr  = str(dtime.hour).zfill(2)
    mi  = str(dtime.minute).zfill(2)
    sec = str(dtime.second).zfill(2)
    msc = str(int(dtime.microsecond/1e4)).zfill(2)
    time = hr + mi + sec + msc
    return date, time


idx=0
for pha in phas:
  if len(pha.split(','))==2:
    # write head line
    ot, lon, lat, _,_,_ = ctlgs[idx].split(',')
    lon = float(lon)
    lon_deg = str(int(lon)).zfill(3)
    lon_min = str(int(100*60*(lon-int(lon)))).zfill(4)
    lat = float(lat)
    lat_deg = str(int(lat)).zfill(2)
    lat_min = str(int(100*60*(lat-int(lat)))).zfill(4)
    ot = UTCDateTime(ot)
    date, time = split_datetime(ot)
    if idx!=0: out.write('\n')
    out.write('{}{} {}{}E{} \n'.format(date+time, lat_deg, lat_min, lon_deg, lon_min))
    idx+=1
  else:
    # write sta line
    net, sta, tp, ts, _,_,_ = pha.split(',')
    tp = UTCDateTime(tp)
    ts = UTCDateTime(ts)
    date = split_datetime(tp)[0]
    ts = ts -(tp - tp.second - tp.microsecond/1e6)
    tp = split_datetime(tp)[1]
    ts = int(100*ts)
    out.write('{:<5}{}  HHZ IPU0{} {}{}{:5}ES 1 \n'.format(sta, net[-2:], date+tp[0:4], tp[4:],' '*7, ts))

