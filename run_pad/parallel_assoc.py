import os, shutil
from obspy import UTCDateTime

# parallel params
pad_dir = '/home/zhouyj/software/PAD'
shutil.copyfile('config.py', os.path.join(pad_dir, 'config.py'))
date_rng = '20150901-20190601'
num_proc = 20
out_root = './output/xls'
ppk_dir = '{}/picks'.format(out_root)
out_pha = '{}/phase_{}.dat'.format(out_root, date_rng)
out_ctlg = '{}/catalog_{}.dat'.format(out_root, date_rng)

# divide by time
ts, te = [UTCDateTime(date) for date in date_rng.split('-')]
dt = (te-ts) / num_proc
for i in range(num_proc):
    t0 = str((ts + i*dt).date)
    t1 = str((ts + (i+1)*dt).date)
    os.system("python {}/run_assoc.py \
        --date_range={} --ppk_dir={} \
        --out_pha={} --out_ctlg={} &"
        .format(pad_dir, date_rng, ppk_dir, out_pha, out_ctlg))
