import os, shutil
from obspy import UTCDateTime

# i/o paths
pal_dir = '/home/zhouyj/software/1_PAL'
shutil.copyfile('config_eg.py', os.path.join(pal_dir, 'config.py'))
sta_file = 'input/example_pal_format1.sta'
data_dir = '/data/Example_data'
out_root = 'output/eg'
out_pick_dir = os.path.join(out_root, 'picks')
# parallel params
time_range = '20190704-20190707'
num_workers = 3

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
    
