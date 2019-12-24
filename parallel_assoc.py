import os, shutil
from obspy import UTCDateTime

# parallel params
pad_dir = '/home/zhouyj/software/PAD'
shutil.copyfile('config.py', os.path.join(pad_dir, 'config.py'))
start_date = '20150901'
end_date   = '20190601'
num_proc = 20
out_root = './output/xls'

# divide by time
ts = UTCDateTime(start_date)
te = UTCDateTime(end_date)
dt = (te-ts) / num_proc
for i in range(num_proc):
    t0 = str((ts + i*dt).date)
    t1 = str((ts + (i+1)*dt).date)
    os.system("python {}/run_assoc.py \
        --time_range={},{} \
        --ppk_dir={}/picks \
        --out_pha={}/phase_{}_{}.dat \
        --out_ctlg={}/catalog_{}_{}.dat &"
        .format(pad_dir, t0,t1, out_root, out_root,t0,t1, out_root))

