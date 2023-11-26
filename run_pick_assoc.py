""" Run picker and associator
    raw waveforms --> picks --> events
"""
import os, glob
import argparse
import numpy as np
from obspy import UTCDateTime
import picker_pal
import associator_pal
import config
import warnings
warnings.filterwarnings("ignore")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str,
                        default='/data/Example_data')
    parser.add_argument('--time_range', type=str,
                        default='20190704-20190707')
    parser.add_argument('--sta_file', type=str,
                        default='input/example_pal_format1.sta')
    parser.add_argument('--out_ctlg', type=str,
                        default='output/eg/catalog.dat')
    parser.add_argument('--out_pha', type=str,
                        default='output/eg/phase.dat')
    parser.add_argument('--out_pick_dir', type=str,
                        default='output/eg/picks')
    args = parser.parse_args()

# PAL config
cfg = config.Config()
get_data_dict = cfg.get_data_dict
read_data = cfg.read_data
sta_dict = cfg.get_sta_dict(args.sta_file)
picker = picker_pal.STA_LTA_Kurtosis(\
    win_sta = cfg.win_sta,
    win_lta = cfg.win_lta,
    trig_thres = cfg.trig_thres,
    p_win = cfg.p_win,
    s_win = cfg.s_win,
    pca_win = cfg.pca_win, 
    pca_range = cfg.pca_range,
    fd_thres = cfg.fd_thres,
    amp_ratio_thres = cfg.amp_ratio_thres,
    amp_win = cfg.amp_win,
    win_kurt = cfg.win_kurt,
    det_gap = cfg.det_gap,
    to_prep = cfg.to_prep,
    freq_band = cfg.freq_band)
associator = associator_pal.PS_Pair_Assoc(\
    sta_dict,
    xy_margin = cfg.xy_margin,
    xy_grid = cfg.xy_grid, 
    z_grids = cfg.z_grids,
    min_sta = cfg.min_sta,
    ot_dev = cfg.ot_dev,
    max_res = cfg.max_res,
    max_drop = cfg.max_drop, 
    vp = cfg.vp)

# i/o paths
out_root = os.path.split(args.out_ctlg)[0]
if not os.path.exists(out_root): os.makedirs(out_root)
if not os.path.exists(args.out_pick_dir): os.makedirs(args.out_pick_dir)
out_ctlg = open(args.out_ctlg,'w')
out_pha = open(args.out_pha,'w')

# get time range
start_date, end_date = [UTCDateTime(date) for date in args.time_range.split('-')]
print('run pick & assoc: raw_waveform --> picks --> events')
print('time range: {} to {}'.format(start_date.date, end_date.date))
# for all days
num_days = (end_date.date - start_date.date).days
for day_idx in range(num_days):
    # get data paths
    date = start_date + day_idx*86400
    data_dict = get_data_dict(date, args.data_dir)
    todel = [net_sta for net_sta in data_dict if net_sta not in sta_dict]
    for net_sta in todel: data_dict.pop(net_sta)
    if data_dict=={}: continue
    # 1. phase picking: waveform --> picks
    fpick_path = os.path.join(args.out_pick_dir, '%s.pick'%date.date)
    out_pick = open(fpick_path,'w')
    for i, st_paths in enumerate(data_dict.values()):
        print('-'*40)
        stream = read_data(st_paths, sta_dict)
        picks_i = picker.pick(stream, out_pick)
        picks = picks_i if i==0 else np.append(picks, picks_i)
    out_pick.close()
    # 2. associate picks: picks --> events
    associator.associate(picks, out_ctlg, out_pha)
out_pha.close()
out_ctlg.close()
