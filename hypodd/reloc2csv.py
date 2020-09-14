""" Format hypoDD output: hyp to csv
"""
from obspy import UTCDateTime
import config

def read_pha(fpha):
    pha_list = []
    f=open(fpha); lines=f.readlines(); f.close()
    for line in lines:
        codes = line.split(',')
        if len(codes)==5:
            ot_name = dtime2str(UTCDateTime(codes[0]))
            pha_list.append([ot_name, []])
        else: pha_list[-1][-1].append(line)
    return pha_list

def dtime2str(dtime):
    date = ''.join(str(dtime).split('T')[0].split('-'))
    time = ''.join(str(dtime).split('T')[1].split(':'))[0:9]
    return date + time

# params
cfg = config.Config()
dep_corr = cfg.dep_corr
out_ctlg = open(cfg.out_ctlg,'w')
out_pha = open(cfg.out_pha,'w')
pha_list = read_pha(cfg.fpha_in)
freloc = 'output/hypoDD.reloc'
f=open(freloc); lines=f.readlines(); f.close()

for line in lines:
    codes = line.split()
    evid = int(codes[0])
    ot_name, pha_lines = pha_list[evid]
    event_name = '{}_{}'.format(evid, ot_name)
    # get loc info
    lat, lon, dep = codes[1:4]
    dep = round(float(dep) - dep_corr, 2)
    mag = codes[16]
    # get time info
    year, mon, day, hour, mnt, sec = codes[10:16]
    sec = '59.999' if sec=='60.000' else sec
    ot = '{}{:0>2}{:0>2}{:0>2}{:0>2}{:0>6}'.format(year, mon, day, hour, mnt, sec)
    out_ctlg.write('{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag))
    out_pha.write('{},{},{},{},{},{}\n'.format(event_name, ot, lat, lon, dep, mag))
    for pha_line in pha_lines: out_pha.write(pha_line)

out_ctlg.close()
out_pha.close()
