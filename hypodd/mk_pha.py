""" Make input phase file for ph2dt
"""
import os
from obspy import UTCDateTime
import config

# i/o paths
cfg = config.Config()
fpha = cfg.fpha_in
fout = open(cfg.fpha_out,'w')
dep_corr = cfg.dep_corr
f=open(fpha); lines=f.readlines(); f.close()

for line in lines:
  codes = line.split(',')
  if len(codes[0])>=14:
    # write head line
    ot = UTCDateTime(codes[0])
    lat, lon, dep, mag = [float(code) for code in codes[1:5]]
    evid = int(codes[-1])
    # format time info
    date = '{:4} {:2} {:2}'.format(ot.year, ot.month, ot.day)
    time = '{:2} {:2} {:5.2f}'.format(ot.hour, ot.minute, ot.second + ot.microsecond/1e6)
    # format loc info
    loc = '{:7.4f} {:9.4f}  {:6.2f} {:4.2f}'.format(lat, lon, dep+dep_corr, mag)
    fout.write('# {} {}  {}  0.00  0.00  0.00  {:>9}\n'.format(date, time, loc, evid))
  else:
    # write sta pick lines
    sta = codes[0].split('.')[1]
    wp, ws = 1., 1.
    if codes[1]!='-1': 
        tp = UTCDateTime(codes[1])
        ttp = tp - ot
        fout.write('{:<5}{}{:6.3f}  {:6.3f}   P\n'.format(sta, ' '*6, ttp, wp))
    if codes[2]!='-1':
        ts = UTCDateTime(codes[2])
        tts = ts - ot
        fout.write('{:<5}{}{:6.3f}  {:6.3f}   S\n'.format(sta, ' '*6, tts, ws))

fout.close()
