""" Make phase file for ph2dt
"""
import os
from obspy import UTCDateTime
import config

# i/o paths
cfg = config.Config()
fpha = cfg.fpha_in
fout = cfg.fpha_out
dep_corr = cfg.dep_corr
f=open(fpha); lines=f.readlines(); f.close()
out=open(fout,'w')

evid = 0
for line in lines:
    codes = line.split(',')
    if len(codes)==5:
        # write head line
        ot, lat, lon, dep, mag = line.split(',')
        # format time info
        ot = UTCDateTime(codes[0])
        date = '{} {:2} {:2}'.format(ot.year, ot.month, ot.day)
        time = '{:2} {:2} {:5.2f}'.format(ot.hour, ot.minute, ot.second + ot.microsecond/1e6)
        # format loc info
        lat, lon, dep, mag = [float(code) for code in codes[1:]]
        loc = '{:7.4f} {:9.4f}  {:6.2f} {:4.2f}'.format(lat, lon, dep+dep_corr, mag)
        out.write('# {} {}  {}  0.00  0.00  0.00  {:>9}\n'.format(date, time, loc, evid))
        evid += 1
    else:
    # write sta pick lines
        sta = codes[1]
        tp, ts = [UTCDateTime(code) for code in codes[2:4]]
        ttp = tp - ot
        tts = ts - ot
        wp, ws = 1., 1.
        out.write('{:<5}{}{:6.3f}  {:6.3f}   P\n'.format(sta, ' '*6, ttp, wp))
        out.write('{:<5}{}{:6.3f}  {:6.3f}   S\n'.format(sta, ' '*6, tts, ws))

out.close()
