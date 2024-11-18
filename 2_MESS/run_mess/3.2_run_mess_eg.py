import os, shutil
from obspy import UTCDateTime

# i/o paths
mess_dir = '/home/zhouyj/software/2_MESS'
data_dir = '/data/Example_data'
time_range = '20190704-20190707'
sta_file = 'input/example_pal_format1.sta'
temp_root = 'output/Example_templates'
temp_pha = 'input/eg_pal.temp'
out_root = 'output/eg'
# divide into time segments
tmin, tmax = [UTCDateTime(ti) for ti in time_range.split('-')]
num_days = 60
num_segs = ((tmax.date - tmin.date).days) // num_days + 1

# run MESS (cpu ver)
shutil.copyfile('config_eg.py', os.path.join(mess_dir, 'config.py'))
for ii in range(num_segs):
  t0, t1 = tmin+ii*num_days*86400, tmin+(ii+1)*num_days*86400
  if t1>tmax: t1 = tmax
  time_range = '-'.join([''.join(str(ti.date).split('-')) for ti in [t0,t1]])
  out_ctlg = os.path.join(out_root, 'catalog_{}.dat'.format(time_range))
  out_pha = os.path.join(out_root, 'phase_{}.dat'.format(time_range))
  os.system("python {}/run_mess.py \
    --data_dir={}  --time_range={} --sta_file={} \
    --temp_root={} --temp_pha={} \
    --out_ctlg={} --out_pha={}"\
    .format(mess_dir, data_dir, time_range, sta_file, temp_root, temp_pha, out_ctlg, out_pha))

