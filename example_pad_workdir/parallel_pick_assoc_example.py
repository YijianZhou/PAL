import os, shutil
from obspy import UTCDateTime

# parallel params
pad_dir = '/home/zhouyj/software/PAD'
shutil.copyfile('config_temp.py', os.path.join(pad_dir, 'config.py'))
data_dir = '/data' # root dir of data
time_range = '20190704-20190717'
sta_file = 'input/example_pad.sta'
num_workers = 13
out_root = 'output'
out_pick_dir = 'output/picks'

# divide by time
start_date, end_date = [UTCDateTime(date) for date in time_range.split('-')]
dt = (end_date - start_date) / num_workers
for proc_idx in range(num_workers):
    t0 = ''.join(str((start_date + proc_idx*dt).date).split('-'))
    t1 = ''.join(str((start_date + (proc_idx+1)*dt).date).split('-'))
    time_range = '{}-{}'.format(t0, t1)
    out_pha = '{}/phase_{}.dat'.format(out_root, time_range)
    out_ctlg = '{}/catalog_{}.dat'.format(out_root, time_range)
    os.system("python {}/run_pick_assoc.py \
        --time_range={} --data_dir={} --sta_file={} \
        --out_pick_dir={} --out_ctlg={} --out_pha={} & " \
        .format(pad_dir, time_range, data_dir, sta_file, out_pick_dir, out_ctlg, out_pha))

