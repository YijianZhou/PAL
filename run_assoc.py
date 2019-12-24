""" Run associator 
    picks --> events
"""
import os, glob
import argparse
import numpy as np
import obspy
from obspy import read, UTCDateTime
import config

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ppk_dir', type=str,
                        default='./output/picks')
    parser.add_argument('--time_range', type=str,
                        default='20180206,20180207')
    parser.add_argument('--out_ctlg', type=str,
                        default='./output/catalog.tmp')
    parser.add_argument('--out_pha', type=str,
                        default='./output/phase.tmp')
    args = parser.parse_args()


# define func
cfg = config.Config()
get_picks = cfg.get_picks
associator = cfg.associator

# i/o paths
out_root = os.path.split(args.out_pha)[0]
if not os.path.exists(out_root): os.makedirs(out_root)
out_ctlg = open(args.out_ctlg,'w')
out_pha = open(args.out_pha,'w')

# get time range
start_date, end_date = [UTCDateTime(date) for date in args.time_range.split(',')]
print('run assoc: picks --> events')
print('time range: {} to {}'.format(start_date, end_date))

# for all days
num_day = (end_date.date - start_date.date).days
for day_idx in range(num_day):

    # get picks
    date = start_date + day_idx*86400
    picks = get_picks(args.ppk_dir, date)

    # 1. associate ot: picks --> events
    event_picks = associator.pick2event(picks)

    # 2. assocuate by locate evnets
    for event_pick in event_picks:
        event_loc, event_pick = associator.locate(event_pick)
        if len(event_loc)==0: continue
        # 4. estimate magnitude
        event_loc_mag = associator.calc_mag(event_pick, event_loc)
        # write catalog and phase
        associator.write_catalog(event_loc_mag, out_ctlg)
        associator.write_phase(event_loc_mag, event_pick, out_pha)

# finish making catalog
out_pha.close()
out_ctlg.close()
