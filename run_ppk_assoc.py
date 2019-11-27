""" Run picker and associator with raw waveforms as input
"""
import os, glob
import argparse
import numpy as np
import obspy
from obspy import read, UTCDateTime
# import PAD package
import config
import data_pipeline as dp
import pickers
import associators

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data3/XJ_SAC/XLS/[1-9]*')
    parser.add_argument('--sta_file', type=str,
                        default='/data3/XJ_SAC/header/station_XLS.dat')
    parser.add_argument('--time_range', type=str,
                        default='20161219,20161222')
    parser.add_argument('--out_ctlg', type=str,
                        default='./output/catalog.tmp')
    parser.add_argument('--out_pha', type=str,
                        default='./output/phase.tmp')
    parser.add_argument('--out_ppk_dir', type=str,
                        default='./output/picks')
    args = parser.parse_args()


# i/o file
out_ctlg = open(args.out_ctlg,'w')
out_pha  = open(args.out_pha, 'w')
sta_dict = dp.get_sta_dict(args.sta_file)
cfg = config.Config()
if not os.path.exists(args.out_ppk_dir): 
    os.makedirs(args.out_ppk_dir)

# define algorithm
picker = pickers.Trad_PS(trig_thres = cfg.trig_thres,
                         s_win = cfg.s_win)
associator = associators.TS_Assoc(sta_dict, cfg.resp_dict,
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
    data_dict = dp.get_xj(args.data_dir, datetime)
    if data_dict=={}: continue
    # set out file
    fpath = os.path.join(args.out_ppk_dir, str(datetime.date)+'.ppk')
    out_ppk = open(fpath,'w')

    # 1. waveform --> phase picks
    # pick all sta
    for i,sta in enumerate(data_dict):
        # read data
        if len(data_dict[sta])<3: continue
        stream  = read(data_dict[sta][0])
        stream += read(data_dict[sta][1])
        stream += read(data_dict[sta][2])

        # phase picking
        picksi = picker.pick(stream, out_ppk)
        if i==0: picks = picksi
        else:    picks = np.append(picks, picksi)

    # 2. associate ot: picks --> events
    event_picks = associator.pick2event(picks)

    # 3. associate by locate evnets
    for event_pick in event_picks:
        event_loc, event_pick = associator.locate(event_pick)
        if len(event_loc)==0: continue
        # 4. estimate magnitude
        event_loc_mag = associator.calc_mag(event_pick, event_loc)
        # write catalog and phase
        associator.write_catalog(event_loc_mag, out_ctlg)
        associator.write_phase(event_loc_mag, event_pick, out_pha)

    out_ppk.close()

# finish making catalog
out_pha.close()
out_ctlg.close()
