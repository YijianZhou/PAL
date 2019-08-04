from obspy import UTCDateTime

net = 'rc'
fpha  = '../output/%s/phase_%s.dat'%(net,net)
fctlg = '../output/%s/catalog_%s.dat'%(net,net)
fout  = 'input/rc.phs'
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
  if len(pha.split(','))==6:
    # write head line
    ot, lon, lat, _,_,_ = ctlgs[idx].split(',')
    lon = abs(float(lon))
    lon_deg = int(lon)
    lon_min = int(100*60*(lon-int(lon)))
    lat = abs(float(lat))
    lat_deg = int(lat)
    lat_min = int(100*60*(lat-int(lat)))
    ot = UTCDateTime(ot)
    date, time = split_datetime(ot)
    if idx!=0: out.write('\n')
    out.write('{}{:0>2} {:0>4}{:0>3}W{:0>4} \n'\
      .format(date+time, lat_deg, lat_min, lon_deg, lon_min))
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
    out.write('{:<5}{}  HHZ IPU0{} {}{}{:5}ES 1 \n'\
      .format(sta, net[-2:], date+tp[0:4], tp[4:],' '*7, ts))

