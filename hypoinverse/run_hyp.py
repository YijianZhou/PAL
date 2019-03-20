import os
import subprocess
fhyp='xj.hyp'
os.system('python mk_sta.py')
os.system('python mk_phs.py')
p = subprocess.Popen(['hyp1.40'], stdin=subprocess.PIPE)
s = "@{}".format(fhyp)
p.communicate(s.encode())
os.system('python sum2csv.py')

