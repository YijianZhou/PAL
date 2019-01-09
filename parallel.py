"""
Simple script for parallel computing
My usage in bash shell:
$screen -S parallel_mkctlg
$source activate obspy
$python parallel.py
then, kill screen to kill the process, if necessary.
"""
import os
from obspy.core import *

# parallel params
start_date = '20160901'
end_date   = '20180901'
num_thread = 10

ts = UTCDateTime(start_date)
te = UTCDateTime(end_date)
dt = (te-ts) /num_thread
for i in range(num_thread):
    t0 = str((ts + i*dt).date)
    t1 = str((ts + (i+1)*dt).date)
    os.system('python mkctlg.py --time_range={},{} \
               --out_pha=./output/phase_{}_{}.dat \
               --out_ctlg=./output/catalog_{}_{}.dat &' \
               .format(t0,t1,t0,t1,t0,t1))

