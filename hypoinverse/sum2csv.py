import glob, os
import numpy as np

def write(out, line, mag):
    info = line.split()
    date = info[0]
    hrmn = info[1]
    sec  = info[2]
    dtime = date+hrmn+sec.zfill(5)
    lat_deg = float(info[3])
    lat_min = float(info[4])
    lat = lat_deg + lat_min/60
    lon_deg = float(line.split('W')[0].split()[-1])
    lon_min = float(line.split('W')[1].split()[0])
    lon = lon_deg + lon_min/60
    dep = float(line.split('W')[1].split()[1])
    mag = float(mags[i].split(',')[-2])
    out.write('{},{:.4f},{:.4f},{:.1f},{}\n'.format(dtime, lat, -lon, dep+2.7, mag))

# i/o paths
net = 'rc'
fnames = glob.glob('output/rc_*.sum')
fmag  = '../output/%s/catalog_%s.dat'%(net,net)
fout1 = 'output/rc.csv'
fout2 = 'output/rc.sum'
fout3 = 'output/rc_bad.csv'
f=open(fmag);  mags=f.readlines();  f.close()
out1 = open(fout1,'w')
out2 = open(fout2,'w')
out3 = open(fout3,'w')

# read sum files
sum_files = []
for fname in fnames:
    f=open(fname); lines=f.readlines(); f.close()
    if len(lines)!=len(mags): continue
    sum_files.append(lines)

for i,mag in enumerate(mags):
  sum_lines = []
  for sum_file in sum_files:
    line = sum_file[i]
    is_loc = 1
    if '-' in line.split() or '#' in line.split(): is_loc = 0
    rms = float(line.split('W')[1].split()[6])
    sum_lines.append((line, is_loc, rms))
  sum_lines = np.array(sum_lines, 
                       dtype=[('line','O'),('is_loc','O'),('rms','O')])
  sum_lines_loc = sum_lines[sum_lines['is_loc']==1]
  num_loc = len(sum_lines_loc)
  if num_loc==0: sum_lines_loc = sum_lines
  sum_lines = np.sort(sum_lines_loc, order='rms')
  write(out1, sum_lines[0]['line'], mag)
  if num_loc==0: write(out3, sum_lines[0]['line'], mag)
  out2.write(sum_lines[0]['line'])
out1.close()
out2.close()
out3.close()
