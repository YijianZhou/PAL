""" Run associator 
    picks --> events
"""
import os, glob
import argparse
import numpy as np
from obspy import UTCDateTime
import config

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ppk_dir', type=str,
                        default='./output/picks')
    parser.add_argument('--date_range', type=str,
                        default='20170224-20170225')
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
out_ctlg = open(args.out_ctlg,'w')
out_pha = open(args.out_pha,'w')

# get date range
start_date, end_date = [UTCDateTime(date) for date in args.date_range.split('-')]
print('run assoc: picks --> events')
print('date range: {} to {}'.format(start_date.date, end_date.date))

# for all days
num_day = (end_date.date - start_date.date).days
for day_idx in range(num_day):
    # 1. get picks
    date = start_date + day_idx*86400
    picks = get_picks(date, args.ppk_dir)
    # 2. associate picks: picks --> event_picks & event_loc
    associator.associate(picks, out_ctlg, out_pha)

# finish making catalog
out_pha.close()
out_ctlg.close()
