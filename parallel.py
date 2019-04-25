import os
from obspy import UTCDateTime

# parallel params
start_date = '20180113' #'20160901'
end_date   = '20180301'
num_thread = 5

# divide by time
ts = UTCDateTime(start_date)
te = UTCDateTime(end_date)
dt = (te-ts) /num_thread
for i in range(num_thread):
    t0 = str((ts + i*dt).date)
    t1 = str((ts + (i+1)*dt).date)
    os.system("python mk_ctlg.py --time_range={},{} \
               --sta_file=/data3/XJ_SAC/header/station_XLS.dat \
               --data_dir=/data3/XJ_SAC/XLS/* \
               --out_pha=./output/xls/phase_{}_{}.dat \
               --out_ctlg=./output/xls/catalog_{}_{}.dat &"\
               .format(t0,t1,t0,t1,t0,t1))

