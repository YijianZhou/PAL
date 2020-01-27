""" Format hypoDD output: hyp to csv
"""
import config

# params
cfg = config.Config()
dep_corr = cfg.dep_corr
out_csv = open(cfg.out_csv,'w')
freloc = 'output/hypoDD.reloc'
f=open(freloc); lines=f.readlines(); f.close()

for line in lines:
    codes = line.split()
    evid = codes[0]
    # get loc info
    lat, lon, dep = codes[1:4]
    dep = round(float(dep) - dep_corr, 2)
    mag = codes[16]
    # get time info
    year, mon, day, hour, mnt, sec = codes[10:16]
    ot = '{}{:0>2}{:0>2}{:0>2}{:0>2}{:0>6}'.format(year, mon, day, hour, mnt, sec)
    out_csv.write('{},{},{},{},{},{}\n'.format(evid, ot, lat, lon, dep, mag))

out_csv.close()
