import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# i/o paths
fctlg = 'output/zsy.csv'; ctlg_name = 'PAD HypoInverse: ZSY Network'
f=open(fctlg); lines=f.readlines(); f.close()
lon_rng = [102.4, 104]
lat_rng = [24.4, 26.6]
dep_rng = [0, 25.]
mag_corr = 0.
dep_corr = 0.
psize = 3.
zoom = 4

lat, lon, dep, mag = [], [], [], []
num=0
for line in lines:
    dtime, lati, loni, depi, mi = line.split(',')
#    if '**********' in [lati, loni]: continue
    lati = float(lati)
    loni = float(loni)
    depi = float(depi) + dep_corr
    if mi=='-inf': continue
    mi = float(mi) + mag_corr
    if loni<lon_rng[0] or loni>lon_rng[1]: continue
    if lati<lat_rng[0] or lati>lat_rng[1]: continue
    if depi<dep_rng[0] or depi>dep_rng[1]: num+=1; continue
    lat.append(lati)
    lon.append(loni)
    dep.append(depi)
    mag.append(mi * psize)
print(num)

# 1. plot depth distribution
fsize=16
plt.figure()
plt.hist(dep)
ax = plt.gca()
plt.xlabel('Depth (km)', fontsize=fsize)
plt.ylabel('Number', fontsize=fsize)
plt.setp(ax.xaxis.get_majorticklabels(), fontsize=fsize)
plt.setp(ax.yaxis.get_majorticklabels(), fontsize=fsize)
plt.title(ctlg_name, fontsize=fsize+2)

# 2. plot loc
plt.style.use('ggplot')
cmap = plt.get_cmap('hot')
xsize = zoom * (lon_rng[1] - lon_rng[0])
ysize = zoom * (lat_rng[1] - lat_rng[0])
fig = plt.figure(figsize=(xsize, ysize))
color = [cmap(1-depi/dep_rng[1]) for depi in dep]
plt.scatter(lon, lat, mag, alpha=0.6, color=color)
# fill up edge
edgex = [lon_rng[0], lon_rng[0], lon_rng[1], lon_rng[1]]
edgey = [lat_rng[0], lat_rng[1], lat_rng[0], lat_rng[1]]
plt.scatter(edgex, edgey, alpha=0)
plt.xlabel('Longitude', fontsize=fsize)
plt.ylabel('Latitude', fontsize=fsize)
plt.title(ctlg_name, fontsize=fsize+2)
# plot colorbar
cb_ax = fig.add_axes([0.76,0.6,0.03,0.25])
cb = mpl.colorbar.ColorbarBase(cb_ax, cmap=cmap)
cb.set_label('Depth (km)')
cticks = np.arange(0,1.1,0.25)
cticklabels = [str((1-ctick)*dep_rng[1]) for ctick in cticks]
cb.set_ticks(cticks)
cb.set_ticklabels(cticklabels)
#plt.savefig(ctlg_path.split('.')[0] + '.pdf')
plt.show()
