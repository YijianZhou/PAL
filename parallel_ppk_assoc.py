import os, shutil
from obspy import UTCDateTime

# parallel params
pad_dir = '/home/zhouyj/software/PAD'
shutil.copyfile('config.py', os.path.join(pad_dir, 'config.py'))
date_rng = ['20150901', '20190601']
num_proc = 20
out_root = './output/xls'

# divide by time
ts = UTCDateTime(date_rng[0])
te = UTCDateTime(date_rng[1])
dt = (te-ts) / num_proc
for i in range(num_proc):
    t0 = str((ts + i*dt).date)
    t1 = str((ts + (i+1)*dt).date)
    os.system("python {}/run_ppk_assoc.py \
        --time_range={},{} \
        --data_dir=/data3/XLS_SAC/*/* \
        --out_pha={}/phase_{}_{}.dat \
        --out_ctlg={}/catalog_{}_{}.dat \
        --out_ppk_dir={}/picks &"
        .format(pad_dir, t0,t1, out_root,t0,t1, out_root,t0,t1, out_root))

