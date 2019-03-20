fpha  = '../output/phase_XJ_XLS.dat'
fctlg = '../output/catalog_XJ_XLS.dat'
fout  = 'xj.phs'
f=open(fpha); phas =f.readlines(); f.close()
f=open(fctlg);ctlgs=f.readlines(); f.close()
out=open(fout,'w')

def split_datetime(datetime):
    date, time = datetime.split('T')
    date = ''.join(date.split('-'))
    t0, t1 = time.split('.')
    time = ''.join(t0.split(':')+[t1[0:2]])
    return date, time

def get_flt(flt):
    return int(100*flt)-100*int(flt)


idx=0
for pha in phas:
    if len(pha.split(','))==2:
        # write head line
        ot, lon, lat, _,_,_ = ctlgs[idx].split(',')
        lon = float(lon)
        lon_int = int(lon)
        lon_flt = 60*(lon-int(lon))
        lon_flt = '{:2}{:2}'.format(int(lon_flt), get_flt(lon_flt))
        lat = float(lat)
        lat_int = int(lat)
        lat_flt = 60*(lat-int(lat))
        lat_flt = '{:2}{:2}'.format(int(lat_flt), get_flt(lat_flt))
        date, time = split_datetime(ot)
        if idx!=0: out.write('\n')
        out.write('{}{} {:4}{}E{:4} \n'.format(date+time, lat_int, lat_flt, lon_int, lon_flt))
        idx+=1
    else:
        # write sta line
        net, sta, tp, ts, _,_,_,_ = pha.split(',')
        date = split_datetime(tp)[0]
        tp = split_datetime(tp)[1]
        ts = split_datetime(ts)[1]
        out.write('{:<5}{}  HHZ IPU0{} {}{} {}ES 0 \n'.format(sta, net[-2:], date+tp[0:4], tp[4:],' '*7, ts[4:]))

