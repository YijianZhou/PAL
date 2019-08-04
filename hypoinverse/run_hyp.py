import os, glob
import numpy as np
import subprocess

# format input
os.system('python mk_sta.py')
os.system('python mk_phs.py')
for fname in glob.glob('output/rc_*.sum'): os.unlink(fname)

for ztri in np.arange(0,40,1):
    # set control file
    fhyp0='rc.hyp'
    fhyp='rc.hyp.tmp'
    f=open(fhyp0); lines=f.readlines(); f.close()
    out=open(fhyp,'w')
    for line in lines:
        info = line.split()
        if len(info)==0: continue
        if info[0] == 'ZTR': line = 'ZTR %s\n'%ztri
        if info[0] == 'SUM': line = "SUM 'output/rc_%s.sum'\n"%ztri
        out.write(line)
    out.close()
    os.rename(fhyp, fhyp0)

    # run hypoinverse
    p = subprocess.Popen(['hyp1.40'], stdin=subprocess.PIPE)
    s = "@{}".format(fhyp0)
    p.communicate(s.encode())

# format output
os.system('python sum2csv.py')
