import os, shutil
from obspy import UTCDateTime

# parallel params
pal_dir = '/home/zhouyj/software/PAL'
shutil.copyfile('config_eg.py', os.path.join(pal_dir, 'config.py'))
data_dir = '/data/Example_data'
time_range = '20190704-20190707'
sta_file = 'input/example_pal.sta'
num_workers = 3
out_root = 'output/eg'
out_pick_dir = 'output/eg/picks'

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
        .format(pal_dir, time_range, data_dir, sta_file, out_pick_dir, out_ctlg, out_pha))
    
