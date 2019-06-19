""" Run detector with picks as input
"""
import os, glob
import argparse
import numpy as np
import obspy
from obspy import read, UTCDateTime
# import PpkDet package
import config
import data_pipeline as dp
import pickers
import detectors

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ppk_dir', type=str,
                        default='./output/picks')
    parser.add_argument('--sta_file', type=str,
                        default='/data3/XJ_SAC/header/station_ZSY.dat')
    parser.add_argument('--time_range', type=str,
                        default='20170914,20170915')
    parser.add_argument('--out_ctlg', type=str,
                        default='./output/catalog.tmp')
    parser.add_argument('--out_pha', type=str,
                        default='./output/phase.tmp')
    args = parser.parse_args()


# i/o file
out_ctlg = open(args.out_ctlg,'w')
out_pha  = open(args.out_pha, 'w')
sta_dict = dp.get_sta_dict(args.sta_file)
cfg = config.Config()

# define algorithm
detector = detectors.TS_Det(sta_dict, cfg.resp_dict,
                            assoc_num = cfg.assoc_num,
                            ot_dev = cfg.ot_dev,
                            ttp_dev = cfg.ttp_dev)

# get time range
start_date = UTCDateTime(args.time_range.split(',')[0])
end_date   = UTCDateTime(args.time_range.split(',')[1])
print('Making catalog')
print('time range: {} to {}'.format(start_date, end_date))

# for all days
num_day = (end_date.date - start_date.date).days
for day_idx in range(num_day):

    # get data paths
    datetime = start_date + day_idx*86400
    picks = dp.get_xj_picks(args.ppk_dir, datetime)

    # 2. associate ot: picks --> events
    event_picks = detector.pick2event(picks)

    # 3. assocuate by locate evnets
    for event_pick in event_picks:
        event_loc, event_pick = detector.locate(event_pick)
        if len(event_loc)==0: continue
        # 4. estimate magnitude
        event_loc_mag = detector.calc_mag(event_pick, event_loc)
        # write catalog and phase
        detector.write_catalog(event_loc_mag, out_ctlg)
        detector.write_phase(event_loc_mag, event_pick, out_pha)

# finish making catalog
out_pha.close()
out_ctlg.close()
