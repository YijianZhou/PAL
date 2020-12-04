""" Format hypoDD output: hyp to csv
"""
from obspy import UTCDateTime
import config

# read fpha with evid
def read_pha(fpha):
    pha_dict = {}
    f=open(fpha); lines=f.readlines(); f.close()
    for line in lines:
        codes = line.split(',')
        if len(codes[0])>=14:
            evid = codes[-1][:-1]
            pha_dict[evid] = []
        else: pha_dict[evid].append(line)
    return pha_dict


# params
cfg = config.Config()
dep_corr = cfg.dep_corr
out_ctlg = open(cfg.out_ctlg,'w')
out_pha = open(cfg.out_pha,'w')
out_pha_all = open(cfg.out_pha_all,'w')
pha_dict = read_pha(cfg.fpha_in)
freloc = 'output/hypoDD.reloc'
f=open(freloc); lines=f.readlines(); f.close()

for line in lines:
    codes = line.split()
    evid = codes[0]
    pha_lines = pha_dict[evid]
    # get loc info
    lat, lon, dep = codes[1:4]
    dep = round(float(dep) - dep_corr, 2)
    mag = codes[16]
    # get time info
    year, mon, day, hour, mnt, sec = codes[10:16]
    sec = '59.999' if sec=='60.000' else sec
    ot = '{}{:0>2}{:0>2}{:0>2}{:0>2}{:0>6}'.format(year, mon, day, hour, mnt, sec)
    out_ctlg.write('{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag))
    out_pha.write('{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag))
    out_pha_all.write('{},{},{},{},{},{}\n'.format(ot, lat, lon, dep, mag, evid))
    for pha_line in pha_lines: 
        out_pha.write(pha_line)
        out_pha_all.write(pha_line)

out_ctlg.close()
out_pha.close()
out_pha_all.close()
