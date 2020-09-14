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

evid = 0
for line in lines:
  codes = line.split(',')
  if len(codes)==5:
    # write head line
    ot, lat, lon, dep, mag = line.split(',')
    # format time info
    ot = UTCDateTime(codes[0])
    date = '{:4} {:2} {:2}'.format(ot.year, ot.month, ot.day)
    time = '{:2} {:2} {:5.2f}'.format(ot.hour, ot.minute, ot.second + ot.microsecond/1e6)
    # format loc info
    lat, lon, dep, mag = [float(code) for code in codes[1:]]
    loc = '{:7.4f} {:9.4f}  {:6.2f} {:4.2f}'.format(lat, lon, dep+dep_corr, mag)
    fout.write('# {} {}  {}  0.00  0.00  0.00  {:>9}\n'.format(date, time, loc, evid))
    evid += 1
  else:
    # write sta pick lines
    sta = codes[0]
    wp, ws = 1., 1.
    if codes[1]!='-1': 
        tp = UTCDateTime(codes[1])
        ttp = tp - ot
        fout.write('{:<5}{}{:6.3f}  {:6.3f}   P\n'.format(sta, ' '*6, ttp, wp))
    if codes[2][:-1]!='-1':
        ts = UTCDateTime(codes[2])
        tts = ts - ot
        fout.write('{:<5}{}{:6.3f}  {:6.3f}   S\n'.format(sta, ' '*6, tts, ws))

fout.close()
