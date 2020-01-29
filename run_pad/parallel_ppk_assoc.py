import os, shutil
from obspy import UTCDateTime

# parallel params
pad_dir = '/home/zhouyj/software/PAD'
shutil.copyfile('config.py', os.path.join(pad_dir, 'config.py'))
data_dir = '/data3/XLS_SAC/*/*'
date_rng = '20150901-20190601'
num_proc = 20
out_root = './output/xls'
out_ppk_dir = '{}/picks'.format(out_root)

# divide by time
start_date, end_date = [UTCDateTime(date) for date in date_rng.split('-')]
dt = (end_date - start_date) / num_proc
for i in range(num_proc):
    t0 = ''.join(str((ts + i*dt).date).split('-'))
    t1 = ''.join(str((ts + (i+1)*dt).date).split('-'))
    date_rng = '{}-{}'.format(t0, t1)
    out_pha = '{}/phase_{}.dat'.format(out_root, date_rng)
    out_ctlg = '{}/catalog_{}.dat'.format(out_root, date_rng)
    os.system("python {}/run_ppk_assoc.py \
        --date_range={} --data_dir={} \
        --out_pha={} --out_ctlg={} --out_ppk_dir={} &"
        .format(pad_dir, date_rng, data_dir, out_pha, out_ctlg, out_ppk_dir))

