""" run HypoInverse (main function)
  Usage:
    (1) modify template hyp control file to tune HypoInverse params
    (2) manually write velo mod (e.g., CRE file), include ref ele if necessary
    (3) set fsta in mk_sta & main
    (4) set net in mk_pha & main, which is the root dir & name code
    (5) python run_hyp.py
  Output:
    csv file: ot, lat, lon, dep, mag
    sum file (hyp)
    arc file (hyp)
"""
import os, glob
import numpy as np
import subprocess

# i/o paths
net='zsy'
ztr_rng = np.arange(0,40,1)
ref_ele = 2.7 # ref ele for CRE mod; be consistant with input .cre file
fhyp_temp = 'temp.hyp'
fhyp = 'input/%s.hyp'%net
f=open(fhyp_temp); lines=f.readlines(); f.close()
fsta = 'input/xj.sta'
fpha = 'input/%s.phs'%net
# if use CRE mod
pmod = 'input/xj_p.cre'
smod = 'input/xj_s.cre'

# format input
os.system('python mk_sta.py')
os.system('python mk_phs.py')
for fname in glob.glob('output/%s_*.sum'%net): os.unlink(fname)

# for all ztr
for ztri in ztr_rng:

    # 1. set control file
    out=open(fhyp,'w')
    for line in lines:
        if line[0:3]=='ZTR': line = "ZTR %s \n"%ztri
        if line[0:3]=='STA': line = "STA '%s' \n"%fsta
        if line[0:3]=='PHS': line = "PHS '%s' \n"%fpha
        if line[0:5]=='CRE 1': line = "CRE 1 '%s' %s T \n"%(pmod, ref_ele)
        if line[0:5]=='CRE 2': line = "CRE 2 '%s' %s T \n"%(smod, ref_ele)
        if line[0:3]=='PRT': line = "PRT 'output/%s_%s.ptr' \n"%(net, ztri)
        if line[0:3]=='ARC': line = "ARC 'output/%s_%s.arc' \n"%(net, ztri)
        if line[0:3]=='SUM': line = "SUM 'output/%s_%s.sum' \n"%(net, ztri)
        out.write(line)
    out.close()

    # 2. run hypoinverse
    p = subprocess.Popen(['hypoInv'], stdin=subprocess.PIPE)
    s = "@{}".format(fhyp)
    p.communicate(s.encode())

# format output
os.system('python sum2csv.py')
for fname in glob.glob('fort.*'): os.unlink(fname)

