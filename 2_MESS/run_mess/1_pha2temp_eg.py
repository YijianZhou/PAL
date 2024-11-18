""" Selection of template events
  Input
    fpha_det: fpha for detection,used to generate evid & event name
    fpha_loc: fpha after location process
  Output
    template phase file
"""
from obspy import UTCDateTime

# i/o paths
fpha_det = 'input/eg_pal.pha'
fpha_loc = 'input/eg_pal_hyp_full.pha'
fout = open('input/eg_pal.temp','w')
# selection criteria
ot_range = '20190704-20190707'
ot_range = [UTCDateTime(code) for code in ot_range.split('-')]
lat_range = [35.5,36.]
lon_range = [-117.8,-117.3]

def dtime2str(dtime):
    date = ''.join(str(dtime).split('T')[0].split('-'))
    time = ''.join(str(dtime).split('T')[1].split(':'))[0:9]
    return date + time

print('get event_dict')
event_dict = {}
evid = 0
f=open(fpha_det); lines=f.readlines(); f.close()
for line in lines:
    codes = line.split(',')
    if len(codes[0])<10: continue
    event_name = dtime2str(UTCDateTime(codes[0]))
    event_dict[str(evid)] = event_name
    evid += 1

print('selecting template')
f=open(fpha_loc); lines=f.readlines(); f.close()
for line in lines:
    codes = line.split(',')
    # if event line
    if len(codes[0])>=14:
        is_temp = True
        ot = UTCDateTime(codes[0])
        if not ot_range[0]<ot<ot_range[1]: is_temp = False
        lat, lon = [float(code) for code in codes[1:3]]
        if not lat_range[0]<lat<lat_range[1]: is_temp = False
        if not lon_range[0]<lon<lon_range[1]: is_temp = False
        evid = codes[-1][:-1]
        event_name = event_dict[evid]
        ot, lat, lon, dep, mag = codes[0:5]
        if is_temp: fout.write('{}_{},{},{},{},{},{}\n'.format(evid, event_name, ot, lat, lon, dep, mag))
    else:
        if is_temp: fout.write(line)
fout.close()
