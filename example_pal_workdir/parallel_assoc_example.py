import os, shutil
from obspy import UTCDateTime

# parallel params
pal_dir = '/home/zhouyj/software/PAL'
shutil.copyfile('config_example.py', os.path.join(pal_dir, 'config.py'))
time_range = '20190704-20190725'
num_workers = 21
out_root = 'output/example'
pick_dir = 'output/example/picks'
sta_file = 'input/example.sta'

# divide by time
start_date, end_date = [UTCDateTime(date) for date in time_range.split('-')]
dt = (end_date - start_date) / num_workers
for proc_idx in range(num_workers):
    t0 = ''.join(str((start_date + proc_idx*dt).date).split('-'))
    t1 = ''.join(str((start_date + (proc_idx+1)*dt).date).split('-'))
    time_range = '{}-{}'.format(t0, t1)
    out_pha = '{}/phase_{}.dat'.format(out_root, time_range)
    out_ctlg = '{}/catalog_{}.dat'.format(out_root, time_range)
    os.system("python {}/run_assoc.py \
        --time_range={} --pick_dir={} --sta_file={} \
        --out_ctlg={} --out_pha={} & " \
        .format(pal_dir, time_range, pick_dir, sta_file, out_ctlg, out_pha))
