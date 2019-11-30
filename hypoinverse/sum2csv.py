""" change HypoInverse output sum file (Hyp71 format) to csv format
"""
import glob, os
import numpy as np

# i/o paths
net = 'zsy'
grd_ele = 2. # typical station elevation
mag_corr = 2. # hypoInv do not support neg mag
lat_code, lon_code = 'N', 'E'
fsums = glob.glob('output/%s_*.sum'%net)
out_csv = open('output/%s.csv'%net,'w')
out_sum = open('output/%s.sum'%net,'w')
out_bad = open('output/%s_bad.csv'%net,'w')
#fmag = '../output/%s/catalog_%s.dat'%(net,net)
#f=open(fmag); mags=f.readlines(); f.close()


def write(out, line):
    codes = line.split()
    date, hrmn, sec = codes[0:3]
    dtime = date + hrmn + sec.zfill(5)
    lat_deg = float(line[20:22])
    lat_min = float(line[23:28])
    lat = lat_deg + lat_min/60
    lon_deg = float(line[29:32])
    lon_min = float(line[33:38])
    lon = lon_deg + lon_min/60
    dep = float(line[38:44])
    mag = float(line[48:52]) - mag_corr
#    mag = float(mags[evid].split(',')[-2])
    out.write('{},{:.4f},{:.4f},{:.1f},{:.1f}\n'.format(dtime, lat, lon, dep+grd_ele, mag))


# read sum files
sum_dict = {}
for fsum in fsums:
  f=open(fsum); lines=f.readlines(); f.close()
  for line in lines:
    evid = line.split()[-1]
    if evid not in sum_dict: sum_dict[evid] = [line]
    else: sum_dict[evid].append(line)


for evid, sum_lines in sum_dict.items():
  # merge sum lines
  sum_list = []
  dtype = [('line','O'),('is_loc','O'),('azm','O'),('npha','O'),('rms','O')]
  for sum_line in sum_lines:
    codes = sum_line.split()
    is_loc = 1 # whether loc reliable
    if '-' in codes or '#' in codes: is_loc = 0
    npha = float(sum_line[52:55])
    azm  = float(sum_line[56:59])
    rms  = float(sum_line[64:69])
    sum_list.append((sum_line, is_loc, azm, npha, rms))
  sum_list = np.array(sum_list, dtype=dtype)
  sum_list_loc = sum_list[sum_list['is_loc']==1]
  num_loc = len(sum_list_loc)
  # if no reliable loc
  if num_loc==0: 
    sum_list_loc = sum_list
    write(out_bad, sum_list[0]['line'])
  # choose best loc
  sum_list_loc = np.sort(sum_list_loc, order=['azm','npha','rms'])
  write(out_csv, sum_list_loc[0]['line'])
  out_sum.write(sum_list_loc[0]['line'])

out_csv.close()
out_sum.close()
out_bad.close()
