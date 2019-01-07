import os, glob
import argparse
import numpy as np
import obspy
from obspy.core import *
# import PpkAssocLoc package
import config
import data_pipeline
import pickers
import associators
import locators

def main(args):

    # i/o file
    if os.path.exists(args.out_ctlg):
        os.unlink(args.out_ctlg)
        os.unlink(args.out_pha)
    out_ctlg = open(args.out_ctlg, 'a')
    out_pha  = open(args.out_pha,  'a')
    sta_dict = data_pipeline.get_sta_dict(args.sta_file)
    cfg = config.Config()

    # define algorithm
    picker     = pickers.Trad_PS()
    associator = associators.Simple_Assoc()
    locator    = locators.Simple_Loc(sta_dict, cfg.resp)

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
        data_dict = data_pipeline.get_zsy(args.data_dir, datetime)

        # 1. waveform --> phase picks
        # pick all sta
        for i,sta in enumerate(data_dict):
            # read data
            stream  = read(data_dict[sta][0])
            stream += read(data_dict[sta][1])
            stream += read(data_dict[sta][2])

            # phase picking
            picksi = picker.pick(stream)
            if i==0: picks = picksi
            else:    picks = np.append(picks, picksi)

        # 2. associate: picks --> events
        event_picks = associator.pick2event(picks)
        # write pahse file
        associator.write(event_picks, out_pha)

        # 3. locate evnets
        for event_pick in event_picks:
            event_loc = locator.locate(event_pick)
            # 4. estimate magnitude
            event_loc_mag = locator.calc_mag(event_pick, event_loc)
            # write catalog
            locator.write(event_loc_mag, out_ctlg)

    # finish making catalog
    out_pha.close()
    out_ctlg.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, 
                        default='/data3/XJ_SAC/ZSY/*')
    parser.add_argument('--sta_file', type=str, 
                        default='/data3/XJ_SAC/header/station.dat')
    parser.add_argument('--time_range', type=str,
                        default='20160901,20180901')
    parser.add_argument('--out_ctlg', type=str,
                        default='./output/catalog.dat')
    parser.add_argument('--out_pha', type=str,
                        default='./output/phase.dat')
    args = parser.parse_args()
    main(args)
